#=
    reach_explicit_blocks!(ϕ, Xhat0, U, n, b, N, overapproximate, blocks, res)

Reachability computation of a given number of two-dimensional blocks of an
affine system with undeterministic inputs.

The variants have the following structure:

INPUT:

- ``ϕ`` -- sparse matrix of a discrete affine system
- ``Xhat0`` -- initial set as a cartesian product over 2d blocks
- ``U`` -- input set of undeterministic inputs
- ``n`` -- ambient dimension
- ``b`` -- number of blocks
- ``N`` -- number of sets computed
- ``overapproximate`` -- function for overapproximation
- ``blocks`` -- the block indices to be computed
- ``res`` -- storage space for the result, a linear array of CartesianProductArray

OUTPUT:

Array of the cartesian product of two-dimensional sets for the given block
indices, and ZeroSet's for the rest of them.
It is obtained by reachability computation of a discrete affine system with
nondeterministic inputs.
=#


# sparse, with input
function reach_explicit_blocks!(ϕ::SparseMatrixCSC{NUM, Int},
                                Xhat0::Vector{<:LazySet{NUM}},
                                U::ConstantNonDeterministicInput,
                                n::Int,
                                b::Int,
                                N::Int,
                                overapproximate::Function,
                                blocks::AbstractVector{Int},
                                res::Vector{CartesianProductArray{NUM}}
                               )::Void where {NUM}
    res[1] = CartesianProductArray(Xhat0)
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int, bj::Int) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]
    @inline G0(bi::Int) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)
    @inline Gk(bi::Int) = ϕpowerk[(2*bi-1):(2*bi), :]

    Xhatk = Vector{LazySet{NUM}}(b)
    Whatk = Vector{LazySet{NUM}}(b)
    dummy_set = ZeroSet(2)
    @inbounds for bi in 1:b
         Xhatk[bi] = dummy_set
    end

    inputs = next_set(U)
    @inbounds for bi in blocks
        Whatk[bi] = overapproximate(G0(bi) * inputs)
    end
    ϕpowerk = copy(ϕ)

    k = 2
    @inbounds while true
        for bi in blocks
            Xhatk_bi = dummy_set
            for bj in 1:b
                if findfirst(F(bi, bj)) != 0
                    Xhatk_bi = Xhatk_bi + F(bi, bj) * Xhat0[bj]
                end
            end
            Xhatk[bi] = overapproximate(Xhatk_bi + Whatk[bi])
        end
        res[k] = CartesianProductArray(copy(Xhatk))

        if k == N
            break
        end

        for bi in blocks
            Whatk[bi] = overapproximate(Whatk[bi] + Gk(bi) * inputs)
        end
        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    nothing
end


# sparse, no input
function reach_explicit_blocks!(ϕ::SparseMatrixCSC{NUM, Int},
                                Xhat0::Vector{<:LazySet{NUM}},
                                n::Int,
                                b::Int,
                                N::Int,
                                overapproximate::Function,
                                blocks::AbstractVector{Int},
                                res::Vector{CartesianProductArray{NUM}}
                               )::Void where {NUM}
    res[1] = CartesianProductArray(Xhat0)
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int, bj::Int) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]

    Xhatk = Vector{LazySet{NUM}}(b)
    dummy_set = ZeroSet(2)
    @inbounds for bi in 1:b
         Xhatk[bi] = dummy_set
    end

    ϕpowerk = copy(ϕ)

    k = 2
    @inbounds while true
        for bi in blocks
            Xhatk_bi = dummy_set
            for bj in 1:b
                if findfirst(F(bi, bj)) != 0
                    Xhatk_bi = Xhatk_bi + F(bi, bj) * Xhat0[bj]
                end
            end
            Xhatk[bi] = overapproximate(Xhatk_bi)
        end
        res[k] = CartesianProductArray(copy(Xhatk))

        if k == N
            break
        end

        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    nothing
end


# dense, with input
function reach_explicit_blocks!(ϕ::AbstractMatrix{NUM},
                                Xhat0::Vector{<:LazySet{NUM}},
                                U::ConstantNonDeterministicInput,
                                n::Int,
                                b::Int,
                                N::Int,
                                overapproximate::Function,
                                blocks::AbstractVector{Int},
                                res::Vector{CartesianProductArray{NUM}}
                               )::Void where {NUM}
    res[1] = CartesianProductArray(Xhat0)
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int, bj::Int) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]
    @inline G0(bi::Int) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)
    @inline Gk(bi::Int) = ϕpowerk[(2*bi-1):(2*bi), :]

    Xhatk = Vector{LazySet{NUM}}(b)
    Whatk = Vector{LazySet{NUM}}(b)
    dummy_set = ZeroSet(2)
    @inbounds for bi in 1:b
         Xhatk[bi] = dummy_set
    end

    inputs = next_set(U)
    @inbounds for bi in blocks
        Whatk[bi] = overapproximate(G0(bi) * inputs)
    end
    ϕpowerk = copy(ϕ)

    k = 2
    arr = Vector{LazySet{NUM}}(b+1)
    @inbounds while true
        for bi in blocks
            for bj in 1:b
                arr[bj] = F(bi, bj) * Xhat0[bj]
            end
            arr[b+1] = Whatk[bi]
            Xhatk[bi] = overapproximate(MinkowskiSumArray(arr))
        end
        res[k] = CartesianProductArray(copy(Xhatk))

        if k == N
            break
        end

        for bi in blocks
            Whatk[bi] = overapproximate(Whatk[bi] + Gk(bi) * inputs)
        end
        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    nothing
end


# dense, no input
function reach_explicit_blocks!(ϕ::AbstractMatrix{NUM},
                                Xhat0::Vector{<:LazySet{NUM}},
                                n::Int,
                                b::Int,
                                N::Int,
                                overapproximate::Function,
                                blocks::AbstractVector{Int},
                                res::Vector{CartesianProductArray{NUM}}
                               )::Void where {NUM}
    res[1] = CartesianProductArray(Xhat0)
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int, bj::Int) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]

    Xhatk = Vector{LazySet{NUM}}(b)
    dummy_set = ZeroSet(2)
    @inbounds for bi in 1:b
         Xhatk[bi] = dummy_set
    end

    ϕpowerk = copy(ϕ)

    k = 2
    arr = Vector{LazySet{NUM}}(b)
    @inbounds while true
        for bi in blocks
            for bj in 1:b
                arr[bj] = F(bi, bj) * Xhat0[bj]
            end
            Xhatk[bi] = overapproximate(MinkowskiSumArray(arr))
        end
        res[k] = CartesianProductArray(copy(Xhatk))

        if k == N
            break
        end

        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    nothing
end


# lazymexp, no input
function reach_explicit_blocks!(ϕ::SparseMatrixExp{NUM},
                                Xhat0::Vector{<:LazySet{NUM}},
                                n::Int,
                                b::Int,
                                N::Int,
                                overapproximate::Function,
                                blocks::AbstractVector{Int},
                                res::Vector{CartesianProductArray{NUM}}
                               )::Void where {NUM}
    res[1] = CartesianProductArray(Xhat0)
    if N == 1
        nothing
        return
    end

    Xhatk = Vector{LazySet{NUM}}(b)
    dummy_set = ZeroSet(2)
    @inbounds for bi in 1:b
         Xhatk[bi] = dummy_set
    end

    ϕpowerk = SparseMatrixExp(ϕ.M)

    k = 2
    @inbounds while true
        for bi in blocks
            ϕpowerk_πbi = get_rows(ϕpowerk, (2*bi-1):(2*bi))
            Xhatk_bi = dummy_set
            for bj in 1:b
                πbi = ϕpowerk_πbi[:, (2*bj-1):(2*bj)]
                if findfirst(πbi) != 0
                    Xhatk_bi = Xhatk_bi + πbi * Xhat0[bj]
                end
            end
            Xhatk[bi] = overapproximate(Xhatk_bi)
        end
        res[k] = CartesianProductArray(copy(Xhatk))

        if k == N
            break
        end

        ϕpowerk.M .= ϕpowerk.M + ϕ.M
        k += 1
    end

    nothing
end


# lazymexp, with input
function reach_explicit_blocks!(ϕ::SparseMatrixExp{NUM},
                                Xhat0::Vector{<:LazySet{NUM}},
                                U::ConstantNonDeterministicInput,
                                n::Int,
                                b::Int,
                                N::Int,
                                overapproximate::Function,
                                blocks::AbstractVector{Int},
                                res::Vector{CartesianProductArray{NUM}}
                               )::Void where {NUM}
    res[1] = CartesianProductArray(Xhat0)
    if N == 1
        nothing
        return
    end

    @inline G0(bi::Int) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)

    Xhatk = Vector{LazySet{NUM}}(b)
    Whatk = Vector{LazySet{NUM}}(b)
    dummy_set = ZeroSet(2)
    @inbounds for bi in 1:b
         Xhatk[bi] = dummy_set
    end

    inputs = next_set(U)
    @inbounds for bi in blocks
        Whatk[bi] = overapproximate(G0(bi) * inputs)
    end
    ϕpowerk = SparseMatrixExp(ϕ.M)

    k = 2
    @inbounds while true
        for bi in blocks
            ϕpowerk_πbi = get_rows(ϕpowerk, (2*bi-1):(2*bi))
            Xhatk_bi = dummy_set
            for bj in 1:b
                πbi = ϕpowerk_πbi[:, (2*bj-1):(2*bj)]
                if findfirst(πbi) != 0
                    Xhatk_bi = Xhatk_bi + πbi * Xhat0[bj]
                end
            end
            Xhatk[bi] = overapproximate(Xhatk_bi + Whatk[bi])
            Whatk[bi] = overapproximate(Whatk[bi] + ϕpowerk_πbi * inputs)
        end
        res[k] = CartesianProductArray(copy(Xhatk))

        if k == N
            break
        end

        ϕpowerk.M .= ϕpowerk.M + ϕ.M
        k += 1
    end

    nothing
end
