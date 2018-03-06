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
                                N::Int,
                                blocks::AbstractVector{Int},
                                partition::AbstractVector{<:AbstractVector{Int}},
                                prop::Property)::Int where {NUM}
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
    end

    @inline G0(bi::AbstractVector{Int}) =
        sparse(1:length(bi), bi, ones(length(bi)), length(bi), n)

    b = length(partition)
    Xhatk = Vector{LazySet{NUM}}(b)
    Whatk = Vector{LazySet{NUM}}(b)
    @inbounds for (i, bi) in enumerate(partition)
         Xhatk[i] = ZeroSet(length(bi))
    end

    inputs = next_set(U)
    @inbounds for i in blocks
        bi = partition[i]
        Whatk[i] = overapproximate(G0(bi) * inputs)
    end
    ϕpowerk = copy(ϕ)

    k = 2
    @inbounds while true
        for i in blocks
            bi = partition[i]
            Xhatk_bi = ZeroSet(length(bi))
            for (j, bj) in enumerate(partition)
                if findfirst(ϕpowerk[bi, bj]) != 0
                    Xhatk_bi = Xhatk_bi + ϕpowerk[bi, bj] * Xhat0[j]
                end
            end
            Xhatk[i] = Xhatk_bi + Whatk[i]
        end
        if !check_property(CartesianProductArray(Xhatk), prop)
            return k
        elseif k == N
            break
        end

        for i in blocks
            bi = partition[i]
            Whatk[i] = overapproximate(Whatk[i] + ϕpowerk[bi, :] * inputs)
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
                                N::Int,
                                blocks::AbstractVector{Int},
                                partition::AbstractVector{<:AbstractVector{Int}},
                                prop::Property)::Int where {NUM}
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
    end

    b = length(partition)
    Xhatk = Vector{LazySet{NUM}}(b)
    @inbounds for (i, bi) in enumerate(partition)
         Xhatk[i] = ZeroSet(length(bi))
    end

    ϕpowerk = copy(ϕ)

    k = 2
    @inbounds while true
        for i in blocks
            bi = partition[i]
            Xhatk_bi = ZeroSet(length(bi))
            for (j, bj) in enumerate(partition)
                if findfirst(ϕpowerk[bi, bj]) != 0
                    Xhatk_bi = Xhatk_bi + ϕpowerk[bi, bj] * Xhat0[j]
                end
            end
            Xhatk[i] = Xhatk_bi
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
                                N::Int,
                                blocks::AbstractVector{Int},
                                partition::AbstractVector{<:AbstractVector{Int}},
                                prop::Property)::Int where {NUM}
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
    end

    @inline G0(bi::AbstractVector{Int}) =
        sparse(1:length(bi), bi, ones(length(bi)), length(bi), n)

    b = length(partition)
    Xhatk = Vector{LazySet{NUM}}(b)
    Whatk = Vector{LazySet{NUM}}(b)
    @inbounds for (i, bi) in enumerate(partition)
         Xhatk[i] = ZeroSet(length(bi))
    end

    inputs = next_set(U)
    @inbounds for i in blocks
        bi = partition[i]
        Whatk[i] = overapproximate(G0(bi) * inputs)
    end
    ϕpowerk = copy(ϕ)

    k = 2
    @inbounds while true
        for i in blocks
            bi = partition[i]
            arr = Vector{LazySet{NUM}}(b+1)
            for (j, bj) in enumerate(partition)
                arr[j] = ϕpowerk[bi, bj] * Xhat0[j]
            end
            arr[b+1] = Whatk[i]
            Xhatk[i] = MinkowskiSumArray(arr)
        end
        if !check_property(CartesianProductArray(Xhatk), prop)
            return k
        elseif k == N
            break
        end

        for i in blocks
            bi = partition[i]
            Whatk[i] = overapproximate(Whatk[i] + ϕpowerk[bi, :] * inputs)
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
                                N::Int,
                                blocks::AbstractVector{Int},
                                partition::AbstractVector{<:AbstractVector{Int}},
                                prop::Property)::Int where {NUM}
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
    end

    b = length(partition)
    Xhatk = Vector{LazySet{NUM}}(b)
    @inbounds for (i, bi) in enumerate(partition)
         Xhatk[i] = ZeroSet(length(bi))
    end

    ϕpowerk = copy(ϕ)

    k = 2
    @inbounds while true
        for i in blocks
            bi = partition[i]
            arr = Vector{LazySet{NUM}}(b+1)
            for (j, bj) in enumerate(partition)
                arr[j] = ϕpowerk[bi, bj] * Xhat0[j]
            end
            Xhatk[i] = MinkowskiSumArray(arr)
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
                                N::Int,
                                blocks::AbstractVector{Int},
                                partition::AbstractVector{<:AbstractVector{Int}},
                                prop::Property)::Int where {NUM}
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
    end

    b = length(partition)
    Xhatk = Vector{LazySet{NUM}}(b)
    @inbounds for (i, bi) in enumerate(partition)
         Xhatk[i] = ZeroSet(length(bi))
    end

    ϕpowerk = SparseMatrixExp(ϕ.M)

    k = 2
    @inbounds while true
        for i in blocks
            bi = partition[i]
            arr = Vector{LazySet{NUM}}(b+1)
            ϕpowerk_πbi = get_rows(ϕpowerk, bi)
            for (j, bj) in enumerate(partition)
                arr[j] = ϕpowerk_πbi[:, bj] * Xhat0[j]
            end
            Xhatk[i] = MinkowskiSumArray(arr)
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
                                N::Int,
                                blocks::AbstractVector{Int},
                                partition::AbstractVector{<:AbstractVector{Int}},
                                prop::Property)::Int where {NUM}
    if !check_property(CartesianProductArray(Xhat0), prop)
        return 1
    elseif N == 1
        return 0
    end

    @inline G0(bi::AbstractVector{Int}) =
        sparse(1:length(bi), bi, ones(length(bi)), length(bi), n)

    b = length(partition)
    Xhatk = Vector{LazySet{NUM}}(b)
    Whatk = Vector{LazySet{NUM}}(b)
    @inbounds for (i, bi) in enumerate(partition)
         Xhatk[i] = ZeroSet(length(bi))
    end

    inputs = next_set(U)
    @inbounds for i in blocks
        bi = partition[i]
        Whatk[i] = overapproximate(G0(bi) * inputs)
    end
    ϕpowerk = SparseMatrixExp(ϕ.M)

    k = 2
    @inbounds while true
        for i in blocks
            bi = partition[i]
            arr = Vector{LazySet{NUM}}(b+1)
            ϕpowerk_πbi = get_rows(ϕpowerk, bi)
            for (j, bj) in enumerate(partition)
                arr[j] = ϕpowerk_πbi[:, bj] * Xhat0[j]
            end
            arr[b+1] = Whatk[i]
            Xhatk[i] = MinkowskiSumArray(arr)
        end
        if !check_property(CartesianProductArray(Xhatk), prop)
            return k
        elseif k == N
            break
        end

        for i in blocks
            ϕpowerk_πbi = get_rows(ϕpowerk, bi)
            Whatk[i] = overapproximate(Whatk[i] + ϕpowerk_πbi * inputs)
        end
        ϕpowerk.M .= ϕpowerk.M + ϕ.M
        k += 1
    end

    return 0
end
