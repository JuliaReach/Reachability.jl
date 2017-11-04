#=
    check_explicit_blocks!(ϕ, Xhat0, U, overapproximate, n, b, N, blocks, prop)

Property checking of a given number of two-dimensional blocks of an affine
system with nondeterministic inputs.

The variants have the following structure:

INPUT:

- ``ϕ`` -- sparse matrix of a discrete affine system
- ``Xhat0`` -- initial set as a cartesian product over 2d blocks
- ``U`` -- input set of undeterministic inputs
- ``overapproximate`` -- function for overapproximation
- ``n`` -- ambient dimension
- ``b`` -- number of blocks
- ``N`` -- number of sets computed
- ``blocks`` -- the block indices to be computed
- ``prop`` -- property to be checked

OUTPUT:

The first time index where the property is violated, and 0 if the property is satisfied.
=#


# sparse, with input
function check_explicit_blocks!(ϕ::SparseMatrixCSC{Float64, Int64},
                                Xhat0::Vector{HPolygon},
                                U::ConstantNonDeterministicInput,
                                overapproximate::Function,
                                n::Int64,
                                b::Int64,
                                N::Int64,
                                blocks::AbstractVector{Int64},
                                prop::Property)::Int64
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
    end

    @inline F(bi::Int64, bj::Int64) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]
    @inline G0(bi::Int64) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)
    @inline Gk(bi::Int64) = ϕpowerk[(2*bi-1):(2*bi), :]

    Xhatk = Vector{LazySet}(b)
    Whatk = Vector{LazySet}(b)
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
            Xhatk[bi] = Xhatk_bi + Whatk[bi]
        end
        if !check_property(CartesianProductArray(Xhatk), prop)
            return k
        elseif k == N
            break
        end

        for bi in blocks
            Whatk[bi] = overapproximate(Whatk[bi] + Gk(bi) * input_state)
        end
        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    return 0
end


# sparse, no input
function check_explicit_blocks!(ϕ::SparseMatrixCSC{Float64, Int64},
                                Xhat0::Vector{HPolygon},
                                n::Int64,
                                b::Int64,
                                N::Int64,
                                blocks::AbstractVector{Int64},
                                prop::Property)::Int64
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
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
            Xhatk[bi] = Xhatk_bi
        end
        if !check_property(CartesianProductArray(Xhatk), prop)
            return k
        elseif k == N
            break
        end

        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    return 0
end


# dense, with input
function check_explicit_blocks!(ϕ::AbstractMatrix{Float64},
                                Xhat0::Vector{HPolygon},
                                U::ConstantNonDeterministicInput,
                                overapproximate::Function,
                                n::Int64,
                                b::Int64,
                                N::Int64,
                                blocks::AbstractVector{Int64},
                                prop::Property)::Int64
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
    end

    @inline F(bi::Int64, bj::Int64) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]
    @inline G0(bi::Int64) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)
    @inline Gk(bi::Int64) = ϕpowerk[(2*bi-1):(2*bi), :]

    Xhatk = Vector{LazySet}(b)
    Whatk = Vector{LazySet}(b)
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
            arr = Vector{LazySet}(b+1)
            for bj in 1:b
                arr[bj] = F(bi, bj) * Xhat0[bj]
            end
            arr[b+1] = Whatk[bi]
            Xhatk[bi] = MinkowskiSumArray(arr)
        end
        if !check_property(CartesianProductArray(Xhatk), prop)
            return k
        elseif k == N
            break
        end

        for bi in blocks
            Whatk[bi] = overapproximate(Whatk[bi] + Gk(bi) * input_state)
        end
        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    return 0
end


# dense, no input
function check_explicit_blocks!(ϕ::AbstractMatrix{Float64},
                                Xhat0::Vector{HPolygon},
                                n::Int64,
                                b::Int64,
                                N::Int64,
                                blocks::AbstractVector{Int64},
                                prop::Property)::Int64
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
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
            arr = Vector{LazySet}(b+1)
            for bj in 1:b
                arr[bj] = F(bi, bj) * Xhat0[bj]
            end
            Xhatk[bi] = MinkowskiSumArray(arr)
        end
        if !check_property(CartesianProductArray(Xhatk), prop)
            return k
        elseif k == N
            break
        end

        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    return 0
end


# lazymexp, no input
function check_explicit_blocks!(ϕ::SparseMatrixExp{Float64},
                                Xhat0::Vector{HPolygon},
                                n::Int64,
                                b::Int64,
                                N::Int64,
                                blocks::AbstractVector{Int64},
                                prop::Property)::Int64
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
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
            arr = Vector{LazySet}(b+1)
            ϕpowerk_πbi = get_rows(ϕpowerk, (2*bi-1):(2*bi))
            for bj in 1:b
                arr[bj] = ϕpowerk_πbi[:, (2*bj-1):(2*bj)] * Xhat0[bj]
            end
            Xhatk[bi] = MinkowskiSumArray(arr)
        end
        if !check_property(CartesianProductArray(Xhatk), prop)
            return k
        elseif k == N
            break
        end

        ϕpowerk.M = ϕpowerk.M + ϕ.M
        k += 1
    end

    return 0
end


# lazymexp, with input
function check_explicit_blocks!(ϕ::SparseMatrixExp{Float64},
                                Xhat0::Vector{HPolygon},
                                U::ConstantNonDeterministicInput,
                                overapproximate::Function,
                                n::Int64,
                                b::Int64,
                                N::Int64,
                                blocks::AbstractVector{Int64},
                                prop::Property)::Int64
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
    end

    @inline G0(bi::Int64) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)

    Xhatk = Vector{LazySet}(b)
    Whatk = Vector{LazySet}(b)
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
            arr = Vector{LazySet}(b+1)
            ϕpowerk_πbi = get_rows(ϕpowerk, (2*bi-1):(2*bi))
            for bj in 1:b
                arr[bj] = ϕpowerk_πbi[:, (2*bj-1):(2*bj)] * Xhat0[bj]
            end
            arr[b+1] = Whatk[bi]
            Xhatk[bi] = MinkowskiSumArray(arr)
        end
        if !check_property(CartesianProductArray(Xhatk), prop)
            return k
        elseif k == N
            break
        end

        for bi in blocks
            ϕpowerk_πbi = get_rows(ϕpowerk, (2*bi-1):(2*bi))
            Whatk[bi] = overapproximate(Whatk[bi] + ϕpowerk_πbi * input_state)
        end
        ϕpowerk.M = ϕpowerk.M + ϕ.M
        k += 1
    end

    return 0
end
