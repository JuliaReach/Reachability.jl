#=
    reach_explicit_block!(ϕ, Xhat0, U, n, b, N, overapproximate, bi, res)

Reachability computation of a given two-dimensional block of an affine system
with undeterministic inputs.

The variants have the following structure:

INPUT:

- ``ϕ`` -- sparse/dense matrix of a discrete affine system
- ``Xhat0`` -- initial set as a cartesian product over 2d blocks
- ``U`` -- input set of undeterministic inputs
- ``n`` -- ambient dimension
- ``b`` -- number of blocks
- ``N`` -- number of sets computed
- ``overapproximate`` -- function for overapproximation
- ``bi`` -- the block index to be computed
- ``res`` -- storage space for the result, a linear array of CartesianProductArray 

OUTPUT:

Array of two-dimensional sets (HPolygons) for the given block index. 
It is obtained by reachability computation of a discrete affine system with
undeterministic inputs, which can be either constant or time-varying.
=#


# sparse, with input
function reach_explicit_block!(ϕ::SparseMatrixCSC{Float64, Int64},
                               Xhat0::Vector{HPolygon},
                               U::ConstantNonDeterministicInput,
                               n::Int64,
                               b::Int64,
                               N::Int64,
                               overapproximate::Function,
                               bi::Int64,
                               res::Vector{HPolygon})
    res[1] = Xhat0[bi]
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int64, bj::Int64) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]
    @inline G0(bi::Int64) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)
    @inline Gk(bi::Int64) = ϕpowerk[(2*bi-1):(2*bi), :]

    inputs = next_set(U)
    Whatk_bi::HPolygon = overapproximate(G0(bi) * inputs)
    ϕpowerk = copy(ϕ)

    voidSet2 = VoidSet(2)

    k = 2
    @inbounds while true
        Xhatk_bi = voidSet2
        for bj in 1:b
            if findfirst(F(bi, bj)) != 0
                Xhatk_bi = Xhatk_bi + F(bi, bj) * Xhat0[bj]
            end
        end
        res[k] = overapproximate(Xhatk_bi + Whatk_bi)

        if k == N
            break
        end

        Whatk_bi = overapproximate(Whatk_bi + Gk(bi) * inputs)
        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    nothing
end


# sparse, no input
function reach_explicit_block!(ϕ::SparseMatrixCSC{Float64, Int64},
                               Xhat0::Vector{HPolygon},
                               n::Int64,
                               b::Int64,
                               N::Int64,
                               overapproximate::Function,
                               bi::Int64,
                               res::Vector{HPolygon})
    res[1] = Xhat0[bi]
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int64, bj::Int64) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]

    ϕpowerk = copy(ϕ)

    voidSet2 = VoidSet(2)

    k = 2
    @inbounds while true
        Xhatk_bi = voidSet2
        for bj in 1:b
            if findfirst(F(bi, bj)) != 0
                Xhatk_bi = Xhatk_bi + F(bi, bj) * Xhat0[bj]
            end
        end
        res[k] = overapproximate(Xhatk_bi)

        if k == N
            break
        end

        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    nothing
end


# dense, with input
function reach_explicit_block!(ϕ::AbstractMatrix{Float64},
                               Xhat0::Vector{HPolygon},
                               U::ConstantNonDeterministicInput,
                               n::Int64,
                               b::Int64,
                               N::Int64,
                               overapproximate::Function,
                               bi::Int64,
                               res::Vector{HPolygon})
    res[1] = Xhat0[bi]
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int64, bj::Int64) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]
    @inline G0(bi::Int64) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)
    @inline Gk(bi::Int64) = ϕpowerk[(2*bi-1):(2*bi), :]

    arr = Vector{LazySet}(b+1)
    inputs = next_set(U)
    arr[b+1] = overapproximate(G0(bi) * inputs)
    ϕpowerk = copy(ϕ)

    k = 2
    @inbounds while true
        for bj in 1:b
            arr[bj] = F(bi, bj) * Xhat0[bj]
        end
        res[k] = overapproximate(MinkowskiSumArray(arr))

        if k == N
            break
        end

        arr[b+1] = overapproximate(arr[b+1] + Gk(bi) * inputs)
        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    nothing
end


# dense, no input
function reach_explicit_block!(ϕ::AbstractMatrix{Float64},
                               Xhat0::Vector{HPolygon},
                               n::Int64,
                               b::Int64,
                               N::Int64,
                               overapproximate::Function,
                               bi::Int64,
                               res::Vector{HPolygon})
    res[1] = Xhat0[bi]
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int64, bj::Int64) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]

    arr = Vector{LazySet}(b)
    ϕpowerk = copy(ϕ)

    k = 2
    @inbounds while true
        for bj in 1:b
            arr[bj] = F(bi, bj) * Xhat0[bj]
        end
        res[k] = overapproximate(MinkowskiSumArray(arr))

        if k == N
            break
        end

        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    nothing
end

# lazymexp, no input
function reach_explicit_block!(ϕ::SparseMatrixExp{Float64},
                               Xhat0::Vector{HPolygon},
                               n::Int64,
                               b::Int64,
                               N::Int64,
                               overapproximate::Function,
                               bi::Int64,
                               res::Vector{HPolygon})
    res[1] = Xhat0[bi]
    if N == 1
        nothing
        return
    end

    voidSet2 = VoidSet(2)
    ϕpowerk = SparseMatrixExp(ϕ.M)

    k = 2
    @inbounds while true
        ϕpowerk_πbi = get_rows(ϕpowerk, (2*bi-1):(2*bi))
        Xhatk_bi = voidSet2
        for bj in 1:b
            πbi = ϕpowerk_πbi[:, (2*bj-1):(2*bj)]
            if findfirst(πbi) != 0
                Xhatk_bi = Xhatk_bi + πbi * Xhat0[bj]
            end
        end
        res[k] = overapproximate(Xhatk_bi)

        if k == N
            break
        end

        ϕpowerk.M = ϕpowerk.M + ϕ.M
        k += 1
    end

    nothing
end

# lazymexp, input
function reach_explicit_block!(ϕ::SparseMatrixExp{Float64},
                               Xhat0::Vector{HPolygon},
                               U::ConstantNonDeterministicInput,
                               n::Int64,
                               b::Int64,
                               N::Int64,
                               overapproximate::Function,
                               bi::Int64,
                               res::Vector{HPolygon})
    res[1] = Xhat0[bi]
    if N == 1
        nothing
        return
    end

    voidSet2 = VoidSet(2)
    inputs = next_set(U)
    Whatk_bi::HPolygon = overapproximate(sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n) * inputs)
    ϕpowerk = SparseMatrixExp(ϕ.M)

    k = 2
    @inbounds while true
        ϕpowerk_πbi = get_rows(ϕpowerk, (2*bi-1):(2*bi))
        Xhatk_bi = voidSet2
        for bj in 1:b
            πbi = ϕpowerk_πbi[:, (2*bj-1):(2*bj)]
            if findfirst(πbi) != 0
                Xhatk_bi = Xhatk_bi + πbi * Xhat0[bj]
            end
        end
        res[k] = overapproximate(Xhatk_bi + Whatk_bi)

        if k == N
            break
        end

        Whatk_bi = overapproximate(Whatk_bi + ϕpowerk_πbi * inputs)
        ϕpowerk.M = ϕpowerk.M + ϕ.M
        k += 1
    end

    nothing
end

