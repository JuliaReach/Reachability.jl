#=
    reach_blocks!(ϕ, Xhat0, U, n, b, termination, overapproximate, blocks, res)

Reachability computation of a given number of two-dimensional blocks of an
affine system with undeterministic inputs.

The variants have the following structure:

### Input

- `ϕ` -- sparse matrix of a discrete affine system
- `Xhat0` -- initial set as a cartesian product over 2d blocks
- `U` -- input set of undeterministic inputs
- `n` -- ambient dimension
- `termination` -- termination check
- `overapproximate` -- function for overapproximation
- `blocks` -- the block indices to be computed
- `partition` -- the partition into blocks
- `res` -- storage space for the result, a linear array of CartesianProductArray

### Output

The index at which the computation has stopped.

### Notes

The reach sets are stored in `res`, an array of the cartesian product for the
given block indices.
=#

# helper functions
@inline proj(bi::UnitRange{Int}, n::Int) =
        sparse(1:length(bi), bi, ones(length(bi)), length(bi), n)
@inline proj(bi::Int, n::Int) = sparse([1], [bi], ones(1), 1, n)
@inline row(ϕpowerk::AbstractMatrix, bi::UnitRange{Int}) = ϕpowerk[bi, :]
@inline row(ϕpowerk::AbstractMatrix, bi::Int) = ϕpowerk[[bi], :]
@inline row(ϕpowerk::SparseMatrixExp, bi::UnitRange{Int}) = get_rows(ϕpowerk, bi)
@inline row(ϕpowerk::SparseMatrixExp, bi::Int) = Matrix(get_row(ϕpowerk, bi))
@inline block(ϕpowerk_πbi::AbstractMatrix, bj::UnitRange{Int}) = ϕpowerk_πbi[:, bj]
@inline block(ϕpowerk_πbi::AbstractMatrix, bj::Int) = ϕpowerk_πbi[:, [bj]]
@inline store!(res, k, X, t0, t1, N) = (res[k] = ReachSet(X, t0, t1))

# sparse
function reach_blocks!(ϕ::SparseMatrixCSC{NUM, Int},
                       Xhat0::Vector{<:LazySet{NUM}},
                       U::Union{ConstantInput, Nothing},
                       overapproximate::Function,
                       overapproximate_inputs::Function,
                       n::Int,
                       N::Int,
                       output_function::Union{Function, Nothing},
                       blocks::AbstractVector{Int},
                       partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
                       δ::NUM,
                       termination::Function,
                       res::Vector{<:ReachSet}
                       )::Tuple{Int, Bool} where {NUM}
    array = CartesianProductArray(Xhat0[blocks])
    X_store = (output_function == nothing) ?
        array :
        box_approximation(output_function(array))
    t0 = zero(δ)
    t1 = δ
    store!(res, 1, X_store, t0, t1, NUM)
    terminate, skip = termination(1, X_store, t0)
    if terminate
        return 1, skip
    end

    b = length(blocks)
    Xhatk = Vector{LazySet{NUM}}(undef, b)
    ϕpowerk = copy(ϕ)

    if U != nothing
        Whatk = Vector{LazySet{NUM}}(undef, b)
        inputs = next_set(U)
        @inbounds for i in 1:b
            bi = partition[blocks[i]]
            Whatk[i] = overapproximate_inputs(1, blocks[i], proj(bi, n) * inputs)
        end
    end

    k = 2
    p = Progress(N, 1, "Computing successors ")
    @inbounds while true
        update!(p, k)
        for i in 1:b
            bi = partition[blocks[i]]
            Xhatk_bi = ZeroSet(length(bi))
            for (j, bj) in enumerate(partition)
                block = ϕpowerk[bi, bj]
                if !iszero(block)
                    Xhatk_bi = Xhatk_bi + block * Xhat0[j]
                end
            end
            Xhatk_bi_lazy = (U == nothing ? Xhatk_bi : Xhatk_bi + Whatk[i])
            Xhatk[i] = (output_function == nothing) ?
                overapproximate(blocks[i], Xhatk_bi_lazy) :
                Xhatk_bi_lazy
        end
        array = CartesianProductArray(copy(Xhatk))
        X_store = (output_function == nothing) ?
            array :
            box_approximation(output_function(array))
        t0 = t1
        t1 += δ
        store!(res, k, X_store, t0, t1, NUM)

        terminate, skip = termination(k, X_store, t0)
        if terminate
            break
        end

        if U != nothing
            for i in 1:b
                bi = partition[blocks[i]]
                Whatk[i] = overapproximate_inputs(k, blocks[i],
                    Whatk[i] + row(ϕpowerk, bi) * inputs)
            end
        end

        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    return k, skip
end

# dense
function reach_blocks!(ϕ::AbstractMatrix{NUM},
                       Xhat0::Vector{<:LazySet{NUM}},
                       U::Union{ConstantInput, Nothing},
                       overapproximate::Function,
                       overapproximate_inputs::Function,
                       n::Int,
                       N::Int,
                       output_function::Union{Function, Nothing},
                       blocks::AbstractVector{Int},
                       partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
                       δ::NUM,
                       termination::Function,
                       res::Vector{<:ReachSet}
                       )::Tuple{Int, Bool} where {NUM}
    array = CartesianProductArray(Xhat0[blocks])
    X_store = (output_function == nothing) ?
        array :
        box_approximation(output_function(array))
    t0 = zero(δ)
    t1 = δ
    store!(res, 1, X_store, t0, t1, NUM)
    terminate, skip = termination(1, X_store, t0)
    if terminate
        return 1, skip
    end

    b = length(blocks)
    Xhatk = Vector{LazySet{NUM}}(undef, b)
    ϕpowerk = copy(ϕ)
    ϕpowerk_cache = similar(ϕ)

    if U != nothing
        Whatk = Vector{LazySet{NUM}}(undef, b)
        inputs = next_set(U)
        @inbounds for i in 1:b
            bi = partition[blocks[i]]
            Whatk[i] = overapproximate_inputs(1, blocks[i], proj(bi, n) * inputs)
        end
    end

    arr_length = (U == nothing) ? length(partition) : length(partition) + 1
    arr = Vector{LazySet{NUM}}(undef, arr_length)
    k = 2
    p = Progress(N, 1, "Computing successors ")
    @inbounds while true
        update!(p, k)
        for i in 1:b
            bi = partition[blocks[i]]
            for (j, bj) in enumerate(partition)
                arr[j] = ϕpowerk[bi, bj] * Xhat0[j]
            end
            if U != nothing
                arr[arr_length] = Whatk[i]
            end
            Xhatk[i] = (output_function == nothing) ?
                overapproximate(blocks[i], MinkowskiSumArray(arr)) :
                MinkowskiSumArray(copy(arr))
        end
        array = CartesianProductArray(copy(Xhatk))

        X_store = (output_function == nothing) ?
            array :
            box_approximation(output_function(array))
        t0 = t1
        t1 += δ
        store!(res, k, X_store, t0, t1, NUM)

        terminate, skip = termination(k, X_store, t0)
        if terminate
            break
        end

        if U != nothing
            for i in 1:b
                bi = partition[blocks[i]]
                Whatk[i] = overapproximate_inputs(k, blocks[i],
                    Whatk[i] + row(ϕpowerk, bi) * inputs)
            end
        end

        _A_mul_B!(ϕpowerk_cache, ϕpowerk, ϕ)
        copyto!(ϕpowerk, ϕpowerk_cache)

        k += 1
    end

    return k, skip
end

# lazy_expm sparse
function reach_blocks!(ϕ::SparseMatrixExp{NUM},
                       assume_sparse::Val{true},
                       Xhat0::Vector{<:LazySet{NUM}},
                       U::Union{ConstantInput, Nothing},
                       overapproximate::Function,
                       overapproximate_inputs::Function,
                       n::Int,
                       N::Int,
                       output_function::Union{Function, Nothing},
                       blocks::AbstractVector{Int},
                       partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
                       δ::NUM,
                       termination::Function,
                       res::Vector{<:ReachSet}
                       )::Tuple{Int, Bool} where {NUM}
    array = CartesianProductArray(Xhat0[blocks])
    X_store = (output_function == nothing) ?
        array :
        box_approximation(output_function(array))
    t0 = zero(δ)
    t1 = δ
    store!(res, 1, X_store, t0, t1, NUM)
    terminate, skip = termination(1, X_store, t0)
    if terminate
        return 1, skip
    end

    b = length(blocks)
    Xhatk = Vector{LazySet{NUM}}(undef, b)
    ϕpowerk = SparseMatrixExp(copy(ϕ.M))

    if U != nothing
        Whatk = Vector{LazySet{NUM}}(undef, b)
        inputs = next_set(U)
        @inbounds for i in 1:b
            bi = partition[blocks[i]]
            Whatk[i] = overapproximate_inputs(1, blocks[i], proj(bi, n) * inputs)
        end
    end

    k = 2
    p = Progress(N, 1, "Computing successors ")
    @inbounds while true
        update!(p, k)
        for i in 1:b
            bi = partition[blocks[i]]
            ϕpowerk_πbi = row(ϕpowerk, bi)
            Xhatk_bi = ZeroSet(length(bi))
            for (j, bj) in enumerate(partition)
                πbi = block(ϕpowerk_πbi, bj)
                if !iszero(πbi)
                    Xhatk_bi = Xhatk_bi + πbi * Xhat0[j]
                end
            end
            Xhatk_bi_lazy = (U == nothing ? Xhatk_bi : Xhatk_bi + Whatk[i])
            Xhatk[i] = (output_function == nothing) ?
                overapproximate(blocks[i], Xhatk_bi_lazy) :
                Xhatk_bi_lazy
            if U != nothing
                Whatk[i] = overapproximate_inputs(k, blocks[i],
                    Whatk[i] + ϕpowerk_πbi * inputs)
            end
        end
        array = CartesianProductArray(copy(Xhatk))
        X_store = (output_function == nothing) ?
            array :
            box_approximation(output_function(array))
        t0 = t1
        t1 += δ
        store!(res, k, X_store, t0, t1, NUM)

        terminate, skip = termination(k, X_store, t0)
        if terminate
            break
        end

        ϕpowerk.M .= ϕpowerk.M + ϕ.M
        k += 1
    end

    return k, skip
end


# lazy_expm dense
function reach_blocks!(ϕ::SparseMatrixExp{NUM},
                       assume_sparse::Val{false},
                       Xhat0::Vector{<:LazySet{NUM}},
                       U::Union{ConstantInput, Nothing},
                       overapproximate::Function,
                       overapproximate_inputs::Function,
                       n::Int,
                       N::Int,
                       output_function::Union{Function, Nothing},
                       blocks::AbstractVector{Int},
                       partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
                       δ::NUM,
                       termination::Function,
                       res::Vector{<:ReachSet}
                       )::Tuple{Int, Bool} where {NUM}
    array = CartesianProductArray(Xhat0[blocks])
    X_store = (output_function == nothing) ?
        array :
        box_approximation(output_function(array))
    t0 = zero(δ)
    t1 = δ
    store!(res, 1, X_store, t0, t1, NUM)
    terminate, skip = termination(1, X_store, t0)
    if terminate
        return 1, skip
    end

    b = length(blocks)
    Xhatk = Vector{LazySet{NUM}}(undef, b)
    ϕpowerk = SparseMatrixExp(copy(ϕ.M))

    if U != nothing
        Whatk = Vector{LazySet{NUM}}(undef, b)
        inputs = next_set(U)
        @inbounds for i in 1:b
            bi = partition[blocks[i]]
            Whatk[i] = overapproximate_inputs(1, blocks[i], proj(bi, n) * inputs)
        end
    end

    arr_length = (U == nothing) ? length(partition) : length(partition) + 1
    arr = Vector{LazySet{NUM}}(undef, arr_length)
    k = 2
    p = Progress(N, 1, "Computing successors ")
    @inbounds while true
        update!(p, k)
        for i in 1:b
            bi = partition[blocks[i]]
            ϕpowerk_πbi = row(ϕpowerk, bi)
            for (j, bj) in enumerate(partition)
                arr[j] = block(ϕpowerk_πbi, bj) * Xhat0[j]
            end
            if U != nothing
                arr[arr_length] = Whatk[i]
            end
            Xhatk[i] = (output_function == nothing) ?
                overapproximate(blocks[i], MinkowskiSumArray(arr)) :
                MinkowskiSumArray(copy(arr))
            if U != nothing
                Whatk[i] = overapproximate_inputs(k, blocks[i],
                    Whatk[i] + ϕpowerk_πbi * inputs)
            end
        end
        array = CartesianProductArray(copy(Xhatk))
        X_store = (output_function == nothing) ?
            array :
            box_approximation(output_function(array))
        t0 = t1
        t1 += δ
        store!(res, k, X_store, t0, t1, NUM)

        terminate, skip = termination(k, X_store, t0)
        if terminate
            break
        end

        ϕpowerk.M .= ϕpowerk.M + ϕ.M
        k += 1
    end

    return k, skip
end
