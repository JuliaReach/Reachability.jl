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

# helper functions
@inline proj(bi::UnitRange{Int}, n::Int) =
        sparse(1:length(bi), bi, ones(length(bi)), length(bi), n)
@inline proj(bi::Int, n::Int) = sparse([1], [bi], ones(1), 1, n)
@inline row(ϕpowerk::AbstractMatrix, bi::UnitRange{Int}) = ϕpowerk[bi, :]
@inline row(ϕpowerk::AbstractMatrix, bi::Int) = ϕpowerk[[bi], :]
@inline row(ϕpowerk::SparseMatrixExp, bi::UnitRange{Int}) =
    get_rows(ϕpowerk, bi)
@inline row(ϕpowerk::SparseMatrixExp, bi::Int) = get_row(ϕpowerk, bi)
@inline block(ϕpowerk_πbi::AbstractMatrix, bj::UnitRange{Int}) =
    ϕpowerk_πbi[:, bj]
@inline block(ϕpowerk_πbi::AbstractMatrix, bj::Int) = ϕpowerk_πbi[:, [bj]]

# sparse
function reach_blocks!(ϕ::SparseMatrixCSC{NUM, Int},
                       Xhat0::Vector{<:LazySet{NUM}},
                       U::Union{ConstantInput, Void},
                       overapproximate::Function,
                       n::Int,
                       N::Int,
                       blocks::AbstractVector{Int},
                       partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
                       res::Vector{CartesianProductArray{NUM}}
                       )::Void where {NUM}
    res[1] = CartesianProductArray(Xhat0[blocks])
    if N == 1
        return nothing
    end

    b = length(blocks)
    Xhatk = Vector{LazySet{NUM}}(b)
    ϕpowerk = copy(ϕ)

    if U != nothing
        Whatk = Vector{LazySet{NUM}}(b)
        inputs = next_set(U)
        @inbounds for i in 1:b
            bi = partition[blocks[i]]
            Whatk[i] = overapproximate(blocks[i], proj(bi, n) * inputs)
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
                if findfirst(block) != 0
                    Xhatk_bi = Xhatk_bi + block * Xhat0[j]
                end
            end
            Xhatk[i] = overapproximate(blocks[i],
                U == nothing ? Xhatk_bi : Xhatk_bi + Whatk[i])
        end
        res[k] = CartesianProductArray(copy(Xhatk))

        if k == N
            break
        end

        if U != nothing
            for i in 1:b
                bi = partition[blocks[i]]
                Whatk[i] = overapproximate(blocks[i],
                                           Whatk[i] + row(ϕpowerk, bi) * inputs)
            end
        end

        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    return nothing
end

# dense
function reach_blocks!(ϕ::AbstractMatrix{NUM},
                       Xhat0::Vector{<:LazySet{NUM}},
                       U::Union{ConstantInput, Void},
                       overapproximate::Function,
                       n::Int,
                       N::Int,
                       blocks::AbstractVector{Int},
                       partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
                       res::Vector{CartesianProductArray{NUM}}
                       )::Void where {NUM}
    res[1] = CartesianProductArray(Xhat0[blocks])
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
            Whatk[i] = overapproximate(blocks[i], proj(bi, n) * inputs)
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
            for (j, bj) in enumerate(partition)
                arr[j] = ϕpowerk[bi, bj] * Xhat0[j]
            end
            if U != nothing
                arr[arr_length] = Whatk[i]
            end
            Xhatk[i] = overapproximate(blocks[i], MinkowskiSumArray(arr))
        end
        res[k] = CartesianProductArray(copy(Xhatk))

        if k == N
            break
        end

        if U != nothing
            for i in 1:b
                bi = partition[blocks[i]]
                Whatk[i] = overapproximate(blocks[i],
                                           Whatk[i] + row(ϕpowerk, bi) * inputs)
            end
        end

        A_mul_B!(ϕpowerk_cache, ϕpowerk, ϕ)
        copy!(ϕpowerk, ϕpowerk_cache)

        k += 1
    end

    return nothing
end

# lazy_expm
function reach_blocks!(ϕ::SparseMatrixExp{NUM},
                       Xhat0::Vector{<:LazySet{NUM}},
                       U::Union{ConstantInput, Void},
                       overapproximate::Function,
                       n::Int,
                       N::Int,
                       blocks::AbstractVector{Int},
                       partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
                       res::Vector{CartesianProductArray{NUM}}
                       )::Void where {NUM}
    res[1] = CartesianProductArray(Xhat0[blocks])
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
            Whatk[i] = overapproximate(blocks[i], proj(bi, n) * inputs)
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
                if findfirst(πbi) != 0
                    Xhatk_bi = Xhatk_bi + πbi * Xhat0[j]
                end
            end
            Xhatk[i] = overapproximate(blocks[i],
                U == nothing ? Xhatk_bi : Xhatk_bi + Whatk[i])
            if U != nothing
                Whatk[i] =
                    overapproximate(blocks[i], Whatk[i] + ϕpowerk_πbi * inputs)
            end
        end
        res[k] = CartesianProductArray(copy(Xhatk))

        if k == N
            break
        end

        ϕpowerk.M .= ϕpowerk.M + ϕ.M
        k += 1
    end

    return nothing
end
