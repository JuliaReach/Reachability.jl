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

Array of the cartesian product of two-dimensional sets (HPolygons) for the
given block indices, and VoidSet's for the rest of them. It is obtained by
reachability computation of a discrete affine system with undeterministic
inputs, which can be either constant or time-varying.
=#


# sparse, with input
function reach_explicit_blocks!(ϕ::SparseMatrixCSC{Float64, Int64},
                                Xhat0::Vector{HPolygon},
                                U::ConstantNonDeterministicInput,
                                n::Int64,
                                b::Int64,
                                N::Int64,
                                overapproximate::Function,
                                blocks::AbstractVector{Int64},
                                res::Vector{CartesianProductArray})
    res[1] = CartesianProductArray(Xhat0)
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int64, bj::Int64) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]
    @inline G0(bi::Int64) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)
    @inline Gk(bi::Int64) = ϕpowerk[(2*bi-1):(2*bi), :]

    Xhatk = Vector{LazySet}(b)
    Whatk = Vector{HPolygon}(b)
    voidSet2 = VoidSet(2)
    @inbounds for bi in 1:b
         Xhatk[bi] = voidSet2
    end

    input_state = start(U).sf
    @inbounds for bi in blocks
        Whatk[bi] = overapproximate(G0(bi) * input_state)
    end
    ϕpowerk = copy(ϕ)

    k = 2
    @inbounds while true
        for bi in blocks
            Xhatk_bi = voidSet2
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
            Whatk[bi] = overapproximate(Whatk[bi] + Gk(bi) * input_state)
        end
        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    nothing
end


# sparse, no input
function reach_explicit_blocks!(ϕ::SparseMatrixCSC{Float64, Int64},
                                Xhat0::Vector{HPolygon},
                                n::Int64,
                                b::Int64,
                                N::Int64,
                                overapproximate::Function,
                                blocks::AbstractVector{Int64},
                                res::Vector{CartesianProductArray})
    res[1] = CartesianProductArray(Xhat0)
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int64, bj::Int64) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]

    Xhatk = Vector{LazySet}(b)
    voidSet2 = VoidSet(2)
    @inbounds for bi in 1:b
         Xhatk[bi] = voidSet2
    end

    ϕpowerk = copy(ϕ)

    k = 2
    @inbounds while true
        for bi in blocks
            Xhatk_bi = voidSet2
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
function reach_explicit_blocks!(ϕ::AbstractMatrix{Float64},
                                Xhat0::Vector{HPolygon},
                                U::ConstantNonDeterministicInput,
                                n::Int64,
                                b::Int64,
                                N::Int64,
                                overapproximate::Function,
                                blocks::AbstractVector{Int64},
                                res::Vector{CartesianProductArray})
    res[1] = CartesianProductArray(Xhat0)
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int64, bj::Int64) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]
    @inline G0(bi::Int64) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)
    @inline Gk(bi::Int64) = ϕpowerk[(2*bi-1):(2*bi), :]

    Xhatk = Vector{LazySet}(b)
    Whatk = Vector{HPolygon}(b)
    voidSet2 = VoidSet(2)
    @inbounds for bi in 1:b
         Xhatk[bi] = voidSet2
    end

    input_state = start(U).sf
    @inbounds for bi in blocks
        Whatk[bi] = overapproximate(G0(bi) * input_state)
    end
    ϕpowerk = copy(ϕ)

    k = 2
    arr = Vector{LazySet}(b+1)
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
            Whatk[bi] = overapproximate(Whatk[bi] + Gk(bi) * input_state)
        end
        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    nothing
end


# dense, no input
function reach_explicit_blocks!(ϕ::AbstractMatrix{Float64},
                                Xhat0::Vector{HPolygon},
                                n::Int64,
                                b::Int64,
                                N::Int64,
                                overapproximate::Function,
                                blocks::AbstractVector{Int64},
                                res::Vector{CartesianProductArray})
    res[1] = CartesianProductArray(Xhat0)
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int64, bj::Int64) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]

    Xhatk = Vector{LazySet}(b)
    voidSet2 = VoidSet(2)
    @inbounds for bi in 1:b
         Xhatk[bi] = voidSet2
    end

    ϕpowerk = copy(ϕ)

    k = 2
    arr = Vector{LazySet}(b)
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
function reach_explicit_blocks!(ϕ::SparseMatrixExp{Float64},
                                Xhat0::Vector{HPolygon},
                                n::Int64,
                                b::Int64,
                                N::Int64,
                                overapproximate::Function,
                                blocks::AbstractVector{Int64},
                                res::Vector{CartesianProductArray})
    res[1] = CartesianProductArray(Xhat0)
    if N == 1
        nothing
        return
    end

    Xhatk = Vector{LazySet}(b)
    voidSet2 = VoidSet(2)
    @inbounds for bi in 1:b
         Xhatk[bi] = voidSet2
    end

    ϕpowerk = SparseMatrixExp(ϕ.M)

    k = 2
    @inbounds while true
        for bi in blocks
            ϕpowerk_πbi = get_rows(ϕpowerk, (2*bi-1):(2*bi))
            Xhatk_bi = voidSet2
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

        ϕpowerk.M = ϕpowerk.M + ϕ.M
        k += 1
    end

    nothing
end


# lazymexp, with input
function reach_explicit_blocks!(ϕ::SparseMatrixExp{Float64},
                                Xhat0::Vector{HPolygon},
                                U::ConstantNonDeterministicInput,
                                n::Int64,
                                b::Int64,
                                N::Int64,
                                overapproximate::Function,
                                blocks::AbstractVector{Int64},
                                res::Vector{CartesianProductArray})
    res[1] = CartesianProductArray(Xhat0)
    if N == 1
        nothing
        return
    end

    @inline G0(bi::Int64) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)

    Xhatk = Vector{LazySet}(b)
    Whatk = Vector{HPolygon}(b)
    voidSet2 = VoidSet(2)
    @inbounds for bi in 1:b
         Xhatk[bi] = voidSet2
    end

    input_state = start(U).sf
    @inbounds for bi in blocks
        Whatk[bi] = overapproximate(G0(bi) * input_state)
    end
    ϕpowerk = SparseMatrixExp(ϕ.M)

    k = 2
    @inbounds while true
        for bi in blocks
            ϕpowerk_πbi = get_rows(ϕpowerk, (2*bi-1):(2*bi))
            Xhatk_bi = voidSet2
            for bj in 1:b
                πbi = ϕpowerk_πbi[:, (2*bj-1):(2*bj)]
                if findfirst(πbi) != 0
                    Xhatk_bi = Xhatk_bi + πbi * Xhat0[bj]
                end
            end
            Xhatk[bi] = overapproximate(Xhatk_bi + Whatk[bi])
            Whatk[bi] = overapproximate(Whatk[bi] + ϕpowerk_πbi * input_state)
        end
        res[k] = CartesianProductArray(copy(Xhatk))

        if k == N
            break
        end

        ϕpowerk.M = ϕpowerk.M + ϕ.M
        k += 1
    end

    nothing
end
