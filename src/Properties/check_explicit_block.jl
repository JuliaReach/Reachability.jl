#=
    check_explicit_block!(ϕ, Xhat0, U, overapproximate, n, b, N, bi, prop)

Property checking of a given two-dimensional block of an affine system with
nondeterministic inputs.

The variants have the following structure:

INPUT:

- ``ϕ`` -- sparse/dense matrix of a discrete affine system
- ``Xhat0`` -- initial set as a cartesian product over 2d blocks
- ``U`` -- input set of undeterministic inputs
- ``overapproximate`` -- function for overapproximation
- ``n`` -- ambient dimension
- ``b`` -- number of blocks
- ``N`` -- number of sets computed
- ``bi`` -- the block index to be computed
- ``prop`` -- property to be checked

OUTPUT:

The first time index where the property is violated, and 0 if the property is satisfied.
=#


# sparse, with input
function check_explicit_block!(ϕ::SparseMatrixCSC{Float64, Int64},
                               Xhat0::Vector{HPolygon},
                               U::ConstantNonDeterministicInput,
                               overapproximate::Function,
                               n::Int64,
                               b::Int64,
                               N::Int64,
                               bi::Int64,
                               prop::Property)::Int64
    if !check_property(Xhat0[bi], prop)
        return 1
    elseif N == 1
        return 0
    end

    @inline F(bi::Int64, bj::Int64) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]
    @inline G0(bi::Int64) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)
    @inline Gk(bi::Int64) = ϕpowerk[(2*bi-1):(2*bi), :]

    inputs = next_set(U)
    Whatk = overapproximate(G0(bi) * inputs)
    ϕpowerk = copy(ϕ)

    dummy_set = ZeroSet(2)

    k = 2
    @inbounds while true
        Xhatk_bi = dummy_set
        for bj in 1:b
            if findfirst(F(bi, bj)) != 0
                Xhatk_bi = Xhatk_bi + F(bi, bj) * Xhat0[bj]
            end
        end
        if !check_property(Xhatk_bi + Whatk, prop)
            return k
        elseif k == N
            break
        end

        Whatk = overapproximate(Whatk + Gk(bi) * inputs)
        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    return 0
end


# sparse, no input
function check_explicit_block!(ϕ::SparseMatrixCSC{Float64, Int64},
                               Xhat0::Vector{HPolygon},
                               n::Int64,
                               b::Int64,
                               N::Int64,
                               bi::Int64,
                               prop::Property)::Int64
    if !check_property(Xhat0[bi], prop)
        return 1
    elseif N == 1
        return 0
    end

    @inline F(bi::Int64, bj::Int64) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]

    dummy_set = ZeroSet(2)

    ϕpowerk = copy(ϕ)

    k = 2
    @inbounds while true
        Xhatk_bi = dummy_set
        for bj in 1:b
            if findfirst(F(bi, bj)) != 0
                Xhatk_bi = Xhatk_bi + F(bi, bj) * Xhat0[bj]
            end
        end
        if !check_property(Xhatk_bi, prop)
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
function check_explicit_block!(ϕ::AbstractMatrix{Float64},
                               Xhat0::Vector{HPolygon},
                               U::ConstantNonDeterministicInput,
                               overapproximate::Function,
                               n::Int64,
                               b::Int64,
                               N::Int64,
                               bi::Int64,
                               prop::Property)::Int64
    if !check_property(Xhat0[bi], prop)
        return 1
    elseif N == 1
        return 0
    end

    @inline F(bi::Int64, bj::Int64) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]
    @inline G0(bi::Int64) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)
    @inline Gk(bi::Int64) = ϕpowerk[(2*bi-1):(2*bi), :]

    arr = Vector{LazySet{Float64}}(b+1)
    inputs = next_set(U)
    arr[b+1] = overapproximate(G0(bi) * inputs)
    ϕpowerk = copy(ϕ)

    k = 2
    @inbounds while true
        for bj in 1:b
            arr[bj] = F(bi, bj) * Xhat0[bj]
        end
        if !check_property(MinkowskiSumArray(arr), prop)
            return k
        elseif k == N
            break
        end

        arr[b+1] = overapproximate(arr[b+1] + Gk(bi) * inputs)
        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    return 0
end


# dense, no input
function check_explicit_block!(ϕ::AbstractMatrix{Float64},
                               Xhat0::Vector{HPolygon},
                               n::Int64,
                               b::Int64,
                               N::Int64,
                               bi::Int64,
                               prop::Property)::Int64
    if !check_property(Xhat0[bi], prop)
        return 1
    elseif N == 1
        return 0
    end

    @inline F(bi::Int64, bj::Int64) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]

    arr = Vector{LazySet{Float64}}(b)
    ϕpowerk = copy(ϕ)

    k = 2
    @inbounds while true
        for bj in 1:b
            arr[bj] = F(bi, bj) * Xhat0[bj]
        end
        if !check_property(MinkowskiSumArray(arr), prop)
            return k
        elseif k == N
            break
        end

        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    return 0
end


# lazymexp, with input
function check_explicit_block!(ϕ::SparseMatrixExp{Float64},
                               Xhat0::Vector{HPolygon},
                               U::ConstantNonDeterministicInput,
                               overapproximate::Function,
                               n::Int64,
                               b::Int64,
                               N::Int64,
                               bi::Int64,
                               prop::Property)::Int64
    if !check_property(Xhat0[bi], prop)
        return 1
    elseif N == 1
        return 0
    end

    @inline G0(bi::Int64) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)

    arr = Vector{LazySet{Float64}}(b+1)
    inputs = next_set(U)
    arr[b+1] = overapproximate(G0(bi) * inputs)

    ϕpowerk = SparseMatrixExp(ϕ.M)

    k = 2
    @inbounds while true
        ϕpowerk_πbi = get_rows(ϕpowerk, (2*bi-1):(2*bi))

        for bj in 1:b
            arr[bj] = ϕpowerk_πbi[:, (2*bj-1):(2*bj)] * Xhat0[bj]
        end
        if !check_property(MinkowskiSumArray(arr), prop)
            return k
        elseif k == N
            break
        end

        arr[b+1] = overapproximate(arr[b+1] + ϕpowerk_πbi * inputs)
        ϕpowerk.M = ϕpowerk.M + ϕ.M
        k += 1
    end

    return 0
end


# lazymexp, no input
function check_explicit_block!(ϕ::SparseMatrixExp{Float64},
                               Xhat0::Vector{HPolygon},
                               n::Int64,
                               b::Int64,
                               N::Int64,
                               bi::Int64,
                               prop::Property)::Int64
    if !check_property(Xhat0[bi], prop)
        return 1
    elseif N == 1
        return 0
    end

    arr = Vector{LazySet{Float64}}(b)

    ϕpowerk = SparseMatrixExp(ϕ.M)

    k = 2
    @inbounds while true
        ϕpowerk_πbi = get_rows(ϕpowerk, (2*bi-1):(2*bi))

        for bj in 1:b
            arr[bj] = ϕpowerk_πbi[:, (2*bj-1):(2*bj)] * Xhat0[bj]
        end
        if !check_property(MinkowskiSumArray(arr), prop)
            return k
        elseif k == N
            break
        end

        ϕpowerk.M = ϕpowerk.M + ϕ.M
        k += 1
    end

    return 0
end


export check_explicit_block
