#=
    reach_blocks!(ϕ, Xhat0, U, n, b, N, overapproximate, blocks, res)

Reachability computation of a given number of two-dimensional blocks of an
affine system with undeterministic inputs.

The variants have the following structure:

INPUT:

- `ϕ` -- sparse matrix of a discrete affine system
- `Xhat0` -- initial set as a cartesian product over 2d blocks
- `U` -- input set of undeterministic inputs
- `n` -- ambient dimension
- `N` -- number of sets computed
- `overapproximate` -- function for overapproximation
- `blocks` -- the block indices to be computed
- `partition` -- the partition into blocks
- `res` -- storage space for the result, a linear array of CartesianProductArray

OUTPUT:

Array of the cartesian product of two-dimensional sets for the given block
indices, and ZeroSet's for the rest of them.
It is obtained by reachability computation of a discrete affine system with
nondeterministic inputs.
=#
#=
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
=#

# sparse
function reach_blocks_parallel!(ϕ::SparseMatrixCSC{NUM, Int},
                       Xhat0::Vector{<:LazySet{NUM}},
                       U::Union{ConstantInput, Void},
                       overapproximate::Function,
                       overapproximate_inputs::Function,
                       n::Int,
                       N::Int,
                       output_function::Union{Function, Void},
                       blocks::AbstractVector{Int},
                       partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
                       res::Vector{OUT}
                       )::Void where {NUM, OUT<:LazySet{NUM}}


   info("Using parallel execution for sparse case")

    array = CartesianProductArray(Xhat0[blocks])
    res[1] = (output_function == nothing) ?
        array :
        box_approximation(output_function(array))
    if N == 1
        return nothing
    end

    b = length(blocks)
    Xhatk = Vector{LazySet{NUM}}(b)
    ϕpowerk = copy(ϕ)

    Whatk = Vector{LazySet{NUM}}(b)

    if U != nothing
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

        Xhatk, Whatk = advection_shared!(k, inputs, ϕpowerk, Xhat0, U, overapproximate, overapproximate_inputs, blocks, output_function, partition, Xhatk, Whatk)

        array = CartesianProductArray(Xhatk)
        res[k] = (output_function == nothing) ?
            array :
            box_approximation(output_function(array))

        if k == N
            break
        end

        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    return nothing
end

function advection_chunk!(
    k::Int,
    inputs,
    ϕpowerk::SparseMatrixCSC{NUM, Int},
    Xhat0::Vector{<:LazySet{NUM}},
    U::Union{ConstantInput, Void},
    overapproximate::Function,
    overapproximate_inputs::Function,
    blocks::AbstractVector{Int},
    output_function::Union{Function, Void},
    partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
    Xhatk::Vector{LazySet{NUM}},
    Whatk::Vector{LazySet{NUM}},
    irange::UnitRange{Int64}) where {NUM}

    XhatkInRange = Vector{LazySet{NUM}}(size(irange))
    WhatkInRange = Vector{LazySet{NUM}}(size(irange))

    count = 1

    for i in irange
        bi = partition[blocks[i]]
        Xhatk_bi = ZeroSet(length(bi))
        for (j, bj) in enumerate(partition)
            block = ϕpowerk[bi, bj]
            if findfirst(block) != 0
                Xhatk_bi = Xhatk_bi + block * Xhat0[j]
            end
        end

        Xhatk_bi_lazy = (U == nothing ? Xhatk_bi : Xhatk_bi + Whatk[i])
        XhatkInRange[count] = ((output_function == nothing) ?
            overapproximate(blocks[i], Xhatk_bi_lazy) :
            Xhatk_bi_lazy)::LazySet{NUM}

        if U != nothing
            bi = partition[blocks[i]]
            WhatkInRange[count] = overapproximate_inputs(k, blocks[i],
                Whatk[i] + row(ϕpowerk, bi) * inputs)
        end
        count += 1
    end

    return XhatkInRange, WhatkInRange

end

# This function retuns the (irange,jrange) indexes assigned to this worker
function myrange(size)
    idx = myid()
    nchunks = length(procs())
    splits = [round(Int, s) for s in linspace(0,size,nchunks+1)]
    splits[idx]+1:splits[idx+1]
end

advection_shared_chunk!(
    k::Int,
    inputs,
    ϕpowerk::SparseMatrixCSC{NUM, Int},
    Xhat0::Vector{<:LazySet{NUM}},
    U::Union{ConstantInput, Void},
    overapproximate::Function,
    overapproximate_inputs::Function,
    blocks::AbstractVector{Int},
    output_function::Union{Function, Void},
    partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
    Xhatk::Vector{LazySet{NUM}},
    Whatk::Vector{LazySet{NUM}}
    ) where {NUM} = advection_chunk!(k, inputs, ϕpowerk, Xhat0, U, overapproximate, overapproximate_inputs, blocks, output_function, partition, Xhatk, Whatk, myrange(length(Xhatk)))

function advection_shared!(
    k::Int,
    inputs,
    ϕpowerk::SparseMatrixCSC{NUM, Int},
    Xhat0::Vector{<:LazySet{NUM}},
    U::Union{ConstantInput, Void},
    overapproximate::Function,
    overapproximate_inputs::Function,
    blocks::AbstractVector{Int},
    output_function::Union{Function, Void},
    partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
    Xhatk::Vector{LazySet{NUM}},
    Whatk::Vector{LazySet{NUM}}) where {NUM}

    tasks = Vector{Future}(length(procs()))

    @sync begin
        for (i, p) in enumerate(procs())
            tasks[i] = remotecall(advection_shared_chunk!, p, k, inputs, ϕpowerk, Xhat0, U, overapproximate, overapproximate_inputs, blocks, output_function, partition, Xhatk, Whatk)
        end
    end

    res = [fetch(t) for t in tasks]
    Xhatk, Whatk = collect(zip(res...))
    return vcat(Xhatk...), vcat(Whatk...)
end


# dense
function  reach_blocks_parallel!(ϕ::AbstractMatrix{NUM},
                       Xhat0::Vector{<:LazySet{NUM}},
                       U::Union{ConstantInput, Void},
                       overapproximate::Function,
                       overapproximate_inputs::Function,
                       n::Int,
                       N::Int,
                       output_function::Union{Function, Void}, # TODO
                       blocks::AbstractVector{Int},
                       partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
                       res::Vector{OUT}
                       )::Void where {NUM, OUT<:LazySet{NUM}}

   info("Using parallel execution for dense case")

    array = CartesianProductArray(Xhat0[blocks])
    res[1] = (output_function == nothing) ?
        array :
        box_approximation(output_function(array))
    if N == 1
        return nothing
    end

    b = length(blocks)
    Xhatk = Vector{LazySet{NUM}}(b)
    ϕpowerk = copy(ϕ)
    ϕpowerk_cache = similar(ϕ)

    if U != nothing
        Whatk = Vector{LazySet{NUM}}(b)
        inputs = next_set(U)
        @inbounds for i in 1:b
            bi = partition[blocks[i]]
            Whatk[i] = overapproximate_inputs(1, blocks[i], proj(bi, n) * inputs)
        end
    end

    arr_length = (U == nothing) ? length(partition) : length(partition) + 1
    arr = Vector{LazySet{NUM}}(arr_length)
    k = 2
    p = Progress(N, 1, "Computing successors ")
    @inbounds while true
        update!(p, k)
        for i in 1:b
            bi = partition[blocks[i]]
            @sync @parallel for (j, bj) in enumerate(partition)
                arr[j] = ϕpowerk[bi, bj] * Xhat0[j]
            end
            if U != nothing
                arr[arr_length] = Whatk[i]

                Whatk[i] = overapproximate_inputs(k, blocks[i],
                    Whatk[i] + row(ϕpowerk, bi) * inputs)
            end

            Xhatk[i] = (output_function == nothing) ?
                overapproximate(blocks[i], MinkowskiSumArray(arr)) :
                MinkowskiSumArray(copy(arr))
        end
        array = CartesianProductArray(copy(Xhatk))
        res[k] = (output_function == nothing) ?
            array :
            box_approximation(output_function(array))

        if k == N
            break
        end

        A_mul_B!(ϕpowerk_cache, ϕpowerk, ϕ)
        copy!(ϕpowerk, ϕpowerk_cache)

        k += 1
    end

    return nothing
end

# lazy_expm sparse
function  reach_blocks_parallel!(ϕ::SparseMatrixExp{NUM},
                       assume_sparse::Val{true},
                       Xhat0::Vector{<:LazySet{NUM}},
                       U::Union{ConstantInput, Void},
                       overapproximate::Function,
                       overapproximate_inputs::Function,
                       n::Int,
                       N::Int,
                       output_function::Union{Function, Void}, # TODO
                       blocks::AbstractVector{Int},
                       partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
                       res::Vector{OUT}
                       )::Void where {NUM, OUT<:LazySet{NUM}}

   info("Using parallel execution for lazy_expm sparse case")
   flush(STDOUT)
    array = CartesianProductArray(Xhat0[blocks])
    res[1] = (output_function == nothing) ?
        array :
        box_approximation(output_function(array))
    if N == 1
        return nothing
    end
    info("AUsing parallel execution for lazy_expm sparse case")
    flush(STDOUT)
    b = length(blocks)
    Xhatk = Vector{LazySet{NUM}}(b)
    ϕpowerk = SparseMatrixExp(copy(ϕ.M))
    info("BUsing parallel execution for lazy_expm sparse case")
    flush(STDOUT)
    if U != nothing
        Whatk = Vector{LazySet{NUM}}(b)
        inputs = next_set(U)
        @inbounds for i in 1:b
            bi = partition[blocks[i]]
            Whatk[i] = overapproximate_inputs(1, blocks[i], proj(bi, n) * inputs)
        end
    end
    info("CUsing parallel execution for lazy_expm sparse case")
    flush(STDOUT)
    k = 2
    p = Progress(N, 1, "Computing successors ")
    @inbounds while true
        update!(p, k)

        Xhatk, Whatk = advection_shared2!(k, inputs, ϕpowerk, Xhat0, U, overapproximate, overapproximate_inputs, blocks, output_function, partition, Xhatk, Whatk)

        array = CartesianProductArray(copy(Xhatk))
        res[k] = (output_function == nothing) ?
            array :
            box_approximation(output_function(array))

        if k == N
            break
        end

        ϕpowerk.M .= ϕpowerk.M + ϕ.M
        k += 1
    end

    return nothing
end

function advection_chunk2!(
    k::Int,
    inputs,
    ϕpowerk::SparseMatrixExp{NUM},
    Xhat0::Vector{<:LazySet{NUM}},
    U::Union{ConstantInput, Void},
    overapproximate::Function,
    overapproximate_inputs::Function,
    blocks::AbstractVector{Int},
    output_function::Union{Function, Void},
    partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
    Xhatk::Vector{LazySet{NUM}},
    Whatk::Vector{LazySet{NUM}},
    irange::UnitRange{Int64}) where {NUM}

    XhatkInRange = Vector{LazySet{NUM}}(size(irange))
    WhatkInRange = Vector{LazySet{NUM}}(size(irange))

    count = 1

    for i in irange
        bi = partition[blocks[i]]
        ϕpowerk_πbi = row(ϕpowerk, bi)
        Xhatk_bi = ZeroSet(length(bi))
        for (j, bj) in enumerate(partition)
            πbi = block(ϕpowerk_πbi, bj)
            if findfirst(πbi) != 0
                Xhatk_bi = Xhatk_bi + πbi * Xhat0[j]
            end
        end
        Xhatk_bi_lazy = (U == nothing ? Xhatk_bi : Xhatk_bi + Whatk[i])
        XhatkInRange[count] = (output_function == nothing) ?
            overapproximate(blocks[i], Xhatk_bi_lazy) :
            Xhatk_bi_lazy
        if U != nothing
            WhatkInRange[count] = overapproximate_inputs(k, blocks[i],
                Whatk[i] + ϕpowerk_πbi * inputs)
        end
        count += 1
    end

    return XhatkInRange, WhatkInRange

end

advection_shared_chunk2!(
    k::Int,
    inputs,
    ϕpowerk::SparseMatrixExp{NUM},
    Xhat0::Vector{<:LazySet{NUM}},
    U::Union{ConstantInput, Void},
    overapproximate::Function,
    overapproximate_inputs::Function,
    blocks::AbstractVector{Int},
    output_function::Union{Function, Void},
    partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
    Xhatk::Vector{LazySet{NUM}},
    Whatk::Vector{LazySet{NUM}}
    ) where {NUM} = advection_chunk2!(k, inputs, ϕpowerk, Xhat0, U, overapproximate, overapproximate_inputs, blocks, output_function, partition, Xhatk, Whatk, myrange(length(Xhatk)))

function advection_shared2!(
    k::Int,
    inputs,
    ϕpowerk::SparseMatrixExp{NUM},
    Xhat0::Vector{<:LazySet{NUM}},
    U::Union{ConstantInput, Void},
    overapproximate::Function,
    overapproximate_inputs::Function,
    blocks::AbstractVector{Int},
    output_function::Union{Function, Void},
    partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
    Xhatk::Vector{LazySet{NUM}},
    Whatk::Vector{LazySet{NUM}}) where {NUM}

    tasks = Vector{Future}(length(procs()))

    @sync begin
        for (i, p) in enumerate(procs())
            tasks[i] = remotecall(advection_shared_chunk2!, p, k, inputs, ϕpowerk, Xhat0, U, overapproximate, overapproximate_inputs, blocks, output_function, partition, Xhatk, Whatk)
        end
    end

    res = [fetch(t) for t in tasks]
    Xhatk, Whatk = collect(zip(res...))
    return vcat(Xhatk...), vcat(Whatk...)
end

# lazy_expm dense
function  reach_blocks_parallel!(ϕ::SparseMatrixExp{NUM},
                       assume_sparse::Val{false},
                       Xhat0::Vector{<:LazySet{NUM}},
                       U::Union{ConstantInput, Void},
                       overapproximate::Function,
                       overapproximate_inputs::Function,
                       n::Int,
                       N::Int,
                       output_function::Union{Function, Void}, # TODO
                       blocks::AbstractVector{Int},
                       partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
                       res::Vector{OUT}
                       )::Void where {NUM, OUT<:LazySet{NUM}}

   info("Using parallel execution for lazy_expm dense case")

    array = CartesianProductArray(Xhat0[blocks])
    res[1] = (output_function == nothing) ?
        array :
        box_approximation(output_function(array))
    if N == 1
        return nothing
    end

    b = length(blocks)
    Xhatk = Vector{LazySet{NUM}}(b)
    ϕpowerk = SparseMatrixExp(copy(ϕ.M))

    if U != nothing
        Whatk = Vector{LazySet{NUM}}(b)
        inputs = next_set(U)
        @inbounds for i in 1:b
            bi = partition[blocks[i]]
            Whatk[i] = overapproximate_inputs(1, blocks[i], proj(bi, n) * inputs)
        end
    end

    arr_length = (U == nothing) ? length(partition) : length(partition) + 1
    arr = Vector{LazySet{NUM}}(arr_length)
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
        res[k] = (output_function == nothing) ?
            array :
            box_approximation(output_function(array))

        if k == N
            break
        end

        ϕpowerk.M .= ϕpowerk.M + ϕ.M
        k += 1
    end

    return nothing
end
