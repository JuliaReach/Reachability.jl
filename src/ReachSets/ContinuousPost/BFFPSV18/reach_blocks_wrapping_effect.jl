# dense
function reach_blocks_wrapping_effect!(
        ϕ::AbstractMatrix{NUM},
        Xhat0::Vector{<:LazySet{NUM}},
        U::Union{ConstantInput, Nothing},
        overapproximate::Function,
        overapproximate_inputs::Function,
        n::Int,
        N::Int,
        output_function::Nothing, # ignored
        blocks::AbstractVector{Int},
        partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
        δ::NUM,
        termination::Function,
        res::Vector{<:AbstractReachSet}
       )::Tuple{Int, Bool} where {NUM}
    X_store = CartesianProductArray(Xhat0)
    t0 = zero(δ)
    t1 = δ
    store!(res, 1, X_store, t0, t1, NUM)
    terminate, skip = termination(1, X_store, t0)
    if terminate
        return 1, skip
    end

    b = length(partition)
    Xhatk = Vector{LazySet{NUM}}(undef, b)

    if U != nothing
        Whatk = Vector{LazySet{NUM}}(undef, b)
        inputs = next_set(U)
        @inbounds for i in 1:b
            bi = partition[blocks[i]]
            Whatk[i] = overapproximate_inputs(1, blocks[i], proj(bi, n) * inputs)
        end
    end

    arr_length = (U == nothing) ? length(partition) : length(partition) + 1
    arr = Vector{LazySet{NUM}}(undef, arr_length)
    k = 2
    p = Progress(N, 1, "Computing successors ")
    @inbounds while true
        ProgressMeter.update!(p, k)
        for i in 1:b
            bi = partition[blocks[i]]
            for (j, bj) in enumerate(partition)
                arr[j] = ϕ[bi, bj] * array(set(res[k-1]))[j]
            end
            if U != nothing
                arr[arr_length] = Whatk[i]
            end
            Xhatk[i] = overapproximate(blocks[i], MinkowskiSumArray(arr))
        end
        X_store = CartesianProductArray(copy(Xhatk))
        t0 = t1
        t1 += δ
        store!(res, k, X_store, t0, t1, NUM)

        terminate, skip = termination(k, X_store, t0)
        if terminate
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

    return k, skip
end
