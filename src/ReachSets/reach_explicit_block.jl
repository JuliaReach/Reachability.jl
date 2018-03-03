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

Array of two-dimensional sets for the given block index. 
It is obtained by reachability computation of a discrete affine system with
nondeterministic inputs.
=#


# sparse, with input
function reach_explicit_block!(ϕ::SparseMatrixCSC{NUM, Int},
                               Xhat0::Vector{<:LazySet{NUM}},
                               U::ConstantNonDeterministicInput,
                               n::Int,
                               b::Int,
                               N::Int,
                               overapproximate::Function,
                               bi::Int,
                               res::Vector{<:LazySet{NUM}}
                              )::Void where {NUM}
    res[1] = Xhat0[bi]
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int, bj::Int) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]
    @inline G0(bi::Int) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)
    @inline Gk(bi::Int) = ϕpowerk[(2*bi-1):(2*bi), :]

    inputs = next_set(U)
    Whatk_bi = overapproximate(G0(bi) * inputs)
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
function reach_explicit_block!(ϕ::SparseMatrixCSC{NUM, Int},
                               Xhat0::Vector{<:LazySet{NUM}},
                               n::Int,
                               b::Int,
                               N::Int,
                               overapproximate::Function,
                               bi::Int,
                               res::Vector{<:LazySet{NUM}}
                              )::Void where {NUM}
    res[1] = Xhat0[bi]
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int, bj::Int) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]

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
function reach_explicit_block!(ϕ::AbstractMatrix{NUM},
                               Xhat0::Vector{<:LazySet{NUM}},
                               U::ConstantNonDeterministicInput,
                               n::Int,
                               b::Int,
                               N::Int,
                               overapproximate::Function,
                               bi::Int,
                               res::Vector{<:LazySet{NUM}}
                              )::Void where {NUM}
    res[1] = Xhat0[bi]
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int, bj::Int) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]
    @inline G0(bi::Int) = sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n)
    @inline Gk(bi::Int) = ϕpowerk[(2*bi-1):(2*bi), :]

    arr = Vector{LazySet{NUM}}(b+1)
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
function reach_explicit_block!(ϕ::AbstractMatrix{NUM},
                               Xhat0::Vector{<:LazySet{NUM}},
                               n::Int,
                               b::Int,
                               N::Int,
                               overapproximate::Function,
                               bi::Int,
                               res::Vector{<:LazySet{NUM}}
                              )::Void where {NUM}
    res[1] = Xhat0[bi]
    if N == 1
        nothing
        return
    end

    @inline F(bi::Int, bj::Int) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]

    arr = Vector{LazySet{NUM}}(b)
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
function reach_explicit_block!(ϕ::SparseMatrixExp{NUM},
                               Xhat0::Vector{<:LazySet{NUM}},
                               n::Int,
                               b::Int,
                               N::Int,
                               overapproximate::Function,
                               bi::Int,
                               res::Vector{<:LazySet{NUM}}
                              )::Void where {NUM}
    res[1] = Xhat0[bi]
    if N == 1
        nothing
        return
    end

    dummy_set = ZeroSet(2)
    ϕpowerk = SparseMatrixExp(ϕ.M)

    k = 2
    @inbounds while true
        ϕpowerk_πbi = get_rows(ϕpowerk, (2*bi-1):(2*bi))
        Xhatk_bi = dummy_set
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
function reach_explicit_block!(ϕ::SparseMatrixExp{NUM},
                               Xhat0::Vector{<:LazySet{NUM}},
                               U::ConstantNonDeterministicInput,
                               n::Int,
                               b::Int,
                               N::Int,
                               overapproximate::Function,
                               bi::Int,
                               res::Vector{<:LazySet{NUM}}
                              )::Void where {NUM}
    res[1] = Xhat0[bi]
    if N == 1
        nothing
        return
    end

    dummy_set = ZeroSet(2)
    inputs = next_set(U)
    Whatk_bi = overapproximate(sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n) * inputs)
    ϕpowerk = SparseMatrixExp(ϕ.M)

    k = 2
    @inbounds while true
        ϕpowerk_πbi = get_rows(ϕpowerk, (2*bi-1):(2*bi))
        Xhatk_bi = dummy_set
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

