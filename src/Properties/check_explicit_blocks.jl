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
function check_explicit_blocks!(ϕ::SparseMatrixCSC{NUM, Int},
                                Xhat0::Vector{<:LazySet{NUM}},
                                U::ConstantNonDeterministicInput,
                                overapproximate::Function,
                                n::Int,
                                b::Int,
                                N::Int,
                                blocks::AbstractVector{Int},
                                partition::AbstractVector{<:AbstractVector{Int}},
                                prop::Property)::Int where {NUM}
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
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
            Xhatk[bi] = Xhatk_bi + Whatk[bi]
        end
        if !check_property(CartesianProductArray(Xhatk), prop)
            return k
        elseif k == N
            break
        end

        for bi in blocks
            Whatk[bi] = overapproximate(Whatk[bi] + Gk(bi) * inputs)
        end
        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    return 0
end


# sparse, no input
function check_explicit_blocks!(ϕ::SparseMatrixCSC{NUM, Int},
                                Xhat0::Vector{<:LazySet{NUM}},
                                n::Int,
                                b::Int,
                                N::Int,
                                blocks::AbstractVector{Int},
                                partition::AbstractVector{<:AbstractVector{Int}},
                                prop::Property)::Int where {NUM}
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
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
function check_explicit_blocks!(ϕ::AbstractMatrix{NUM},
                                Xhat0::Vector{<:LazySet{NUM}},
                                U::ConstantNonDeterministicInput,
                                overapproximate::Function,
                                n::Int,
                                b::Int,
                                N::Int,
                                blocks::AbstractVector{Int},
                                partition::AbstractVector{<:AbstractVector{Int}},
                                prop::Property)::Int where {NUM}
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
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
            arr = Vector{LazySet{NUM}}(b+1)
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
            Whatk[bi] = overapproximate(Whatk[bi] + Gk(bi) * inputs)
        end
        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    return 0
end


# dense, no input
function check_explicit_blocks!(ϕ::AbstractMatrix{NUM},
                                Xhat0::Vector{<:LazySet{NUM}},
                                n::Int,
                                b::Int,
                                N::Int,
                                blocks::AbstractVector{Int},
                                partition::AbstractVector{<:AbstractVector{Int}},
                                prop::Property)::Int where {NUM}
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
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
            arr = Vector{LazySet{NUM}}(b+1)
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
function check_explicit_blocks!(ϕ::SparseMatrixExp{NUM},
                                Xhat0::Vector{<:LazySet{NUM}},
                                n::Int,
                                b::Int,
                                N::Int,
                                blocks::AbstractVector{Int},
                                partition::AbstractVector{<:AbstractVector{Int}},
                                prop::Property)::Int where {NUM}
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
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
            arr = Vector{LazySet{NUM}}(b+1)
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

        ϕpowerk.M .= ϕpowerk.M + ϕ.M
        k += 1
    end

    return 0
end


# lazymexp, with input
function check_explicit_blocks!(ϕ::SparseMatrixExp{NUM},
                                Xhat0::Vector{<:LazySet{NUM}},
                                U::ConstantNonDeterministicInput,
                                overapproximate::Function,
                                n::Int,
                                b::Int,
                                N::Int,
                                blocks::AbstractVector{Int},
                                partition::AbstractVector{<:AbstractVector{Int}},
                                prop::Property)::Int where {NUM}
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
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
            arr = Vector{LazySet{NUM}}(b+1)
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
            Whatk[bi] = overapproximate(Whatk[bi] + ϕpowerk_πbi * inputs)
        end
        ϕpowerk.M .= ϕpowerk.M + ϕ.M
        k += 1
    end

    return 0
end
