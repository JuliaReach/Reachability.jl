# dense
function reach_blocks_wrapping_effect!(
        ϕ::AbstractMatrix{NUM},
        Xhat0::Vector{<:LazySet{NUM}},
        U::Union{ConstantInput, Void},
        overapproximate::Function,
        overapproximate_inputs::Function,
        n::Int,
        N::Int,
        output_function::Void, # ignored
        blocks::AbstractVector{Int},
        partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
        res::Vector{OUT}
       )::Void where {NUM, OUT<:LazySet{NUM}}
    res[1] = CartesianProductArray(Xhat0)
    if N == 1
        return nothing
    end

    b = length(partition)
    Xhatk = Vector{LazySet{NUM}}(b)

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
            for (j, bj) in enumerate(partition)
                arr[j] = ϕ[bi, bj] * array(res[k-1])[j]
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
                Whatk[i] =
                    overapproximate_inputs(k, blocks[i], row(ϕ, bi) * inputs)
            end
        end

        k += 1
    end

    return nothing
end
