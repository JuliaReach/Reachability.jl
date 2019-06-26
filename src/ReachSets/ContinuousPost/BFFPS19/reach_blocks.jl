

# sparse
function reach_blocks!(ϕ::SparseMatrixCSC{NUM, Int},
                       Xhat0::Vector{<:LazySet{NUM}},
                       U::Union{ConstantInput, Nothing},
                       overapproximate::Function,
                       overapproximate_inputs::Function,
                       n::Int,
                       N::Union{Int, Nothing},
                       output_function::Union{Function, Nothing},
                       blocks::AbstractVector{Int},
                       partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
                       δ::NUM,
                       guards_proj::Vector{<:LazySet{NUM}},
                       block_options,
                       vars::Vector{Int},
                       termination::Function,
                       progress_meter::Union{Progress, Nothing},
                       res::Vector{<:ReachSet}
                       )::Tuple{Int, Bool} where {NUM}
    update!(progress_meter, 1)
    array = CartesianProductArray(Xhat0[blocks])
    X_store = (output_function == nothing) ?
        array :
        box_approximation(output_function(array))
    t0 = zero(δ)
    t1 = δ

    diff_blocks = setdiff(collect(1:length(partition)),blocks)

    b = length(blocks)
    bd = length(diff_blocks)

    array_d = CartesianProductArray(Xhat0[diff_blocks])
    X_store_d = (output_function == nothing) ?
            array_d :
            box_approximation(output_function(array_d))
    X_store = combine_cpas(X_store, X_store_d, blocks, diff_blocks)
    store!(res, 1, X_store, t0, t1, N)
    terminate, skip = termination(1, X_store, t0)
    if terminate
        return 1, skip
    end

    b = length(blocks)
    Xhatk = Vector{LazySet{NUM}}(undef, b)
    Xhatk_d = Vector{LazySet{NUM}}(undef, bd)
    ϕpowerk = copy(ϕ)

    if U != nothing
        Whatk_blocks = Vector{LazySet{NUM}}(undef, b)
        inputs = next_set(U)
        @inbounds for i in 1:b
            bi = partition[blocks[i]]
            Whatk_blocks[i] = overapproximate_inputs(1, blocks[i], proj(bi, n) * inputs)
       end
    end

    if U != nothing
       Whatk_diff_blocks = Vector{LazySet{NUM}}(undef, bd)
       inputs = next_set(U)
       @inbounds for i in 1:bd
           bi = partition[diff_blocks[i]]
           Whatk_diff_blocks[i] = overapproximate_inputs(1, diff_blocks[i], proj(bi, n) * inputs)
        end
    end

    k = 2
    @inbounds while true
        update!(progress_meter, k)
        for i in 1:b
            bi = partition[blocks[i]]
            Xhatk_bi = ZeroSet(length(bi))
            for (j, bj) in enumerate(partition)
                block = ϕpowerk[bi, bj]
                if !iszero(block)
                    Xhatk_bi = Xhatk_bi + block * Xhat0[j]
                end
            end
            Xhatk_bi_lazy = (U == nothing ? Xhatk_bi : Xhatk_bi + Whatk_blocks[i])
            Xhatk[i] = (output_function == nothing) ?
                overapproximate(blocks[i], Xhatk_bi_lazy) :
                Xhatk_bi_lazy
        end
        array = CartesianProductArray(copy(Xhatk))
        X_store = (output_function == nothing) ?
            array :
            box_approximation(output_function(array))


        terminate, skip, rs = termination(k, X_store, t0)

        if  !(isdisjoint(X_store, UnionSetArray(guards_proj)))
            for i in 1:bd
                bi = partition[diff_blocks[i]]
                Xhatk_bi = ZeroSet(length(bi))
                for (j, bj) in enumerate(partition)
                    block = ϕpowerk[bi, bj]
                    if !iszero(block)
                        Xhatk_bi = Xhatk_bi + block * Xhat0[j]
                    end
                end
                Xhatk_bi_lazy = (U == nothing ? Xhatk_bi : Xhatk_bi + Whatk_diff_blocks[i])
                Xhatk[i] = (output_function == nothing) ?
                    overapproximate(diff_blocks[i], Xhatk_bi_lazy) :
                    Xhatk_bi_lazy
            end

             array_d = CartesianProductArray(copy(Xhatk_d))

             X_store_d = (output_function == nothing) ?
                            array_d :
                            box_approximation(output_function(array_d))
            if !terminate
                X_store = combine_cpas(LazySets.Approximations.overapproximate(rs, CartesianProductArray, block_options), X_store_d, blocks, diff_blocks)
            end
        end

        t0 = t1
        t1 += δ
        store!(res, k, X_store, t0, t1, N)

        terminate, skip = termination(k, X_store, t0)
        if terminate
            break
        end

        if U != nothing
            for i in 1:b
                bi = partition[blocks[i]]
                Whatk_blocks[i] = overapproximate_inputs(k, blocks[i],
                    Whatk_blocks[i] + row(ϕpowerk, bi) * inputs)
            end
        end

         ##for diffs
        if U != nothing
            for i in 1:bd
                bi = partition[diff_blocks[i]]
                Whatk_diff_blocks[i] = overapproximate_inputs(k, diff_blocks[i],
                    Whatk_diff_blocks[i] + row(ϕpowerk, bi) * inputs)
            end
        end

        ϕpowerk = ϕpowerk * ϕ
        k += 1
    end

    return k, skip
end

# dense
function reach_blocks!(ϕ::AbstractMatrix{NUM},
                       Xhat0::Vector{<:LazySet{NUM}},
                       U::Union{ConstantInput, Nothing},
                       overapproximate::Function,
                       overapproximate_inputs::Function,
                       n::Int,
                       N::Union{Int, Nothing},
                       output_function::Union{Function, Nothing},
                       blocks::AbstractVector{Int},
                       partition::AbstractVector{<:Union{AbstractVector{Int}, Int}},
                       δ::NUM,
                       guards_proj::Vector{<:LazySet{NUM}},
                       block_options,
                       vars::Vector{Int},
                       termination::Function,
                       progress_meter::Union{Progress, Nothing},
                       res::Vector{<:ReachSet}
                       )::Tuple{Int, Bool} where {NUM}
    update!(progress_meter, 1)
    array = CartesianProductArray(Xhat0[blocks])
    X_store = (output_function == nothing) ?
        array :
        box_approximation(output_function(array))
    t0 = zero(δ)
    t1 = δ

    diff_blocks = setdiff(collect(1:length(partition)),blocks)

    b = length(blocks)
    bd = length(diff_blocks)

    array_d = CartesianProductArray(Xhat0[diff_blocks])
    X_store_d = (output_function == nothing) ?
        array_d :
        box_approximation(output_function(array_d))

    terminate, skip = termination(1, X_store, t0)
    X_store = combine_cpas(X_store, X_store_d, blocks, diff_blocks)
    store!(res, 1, X_store, t0, t1, N)
    if terminate
        return 1, skip
    end

    Xhatk = Vector{LazySet{NUM}}(undef, b)
    Xhatk_d = Vector{LazySet{NUM}}(undef, bd)
    ϕpowerk = copy(ϕ)
    ϕpowerk_cache = similar(ϕ)

    if U != nothing
        Whatk_blocks = Vector{LazySet{NUM}}(undef, b)
        inputs = next_set(U)
        @inbounds for i in 1:b
            bi = partition[blocks[i]]
            Whatk_blocks[i] = overapproximate_inputs(1, blocks[i], proj(bi, n) * inputs)
        end
    end

     if U != nothing
        Whatk_diff_blocks = Vector{LazySet{NUM}}(undef, bd)
        inputs = next_set(U)
        @inbounds for i in 1:bd
            bi = partition[diff_blocks[i]]
            Whatk_diff_blocks[i] = overapproximate_inputs(1, diff_blocks[i], proj(bi, n) * inputs)
        end
    end

    arr_length = (U == nothing) ? length(partition) : length(partition) + 1
    arr = Vector{LazySet{NUM}}(undef, arr_length)
    k = 2
    @inbounds while true
        update!(progress_meter, k)
        for i in 1:b
            bi = partition[blocks[i]]
            for (j, bj) in enumerate(partition)
                arr[j] = ϕpowerk[bi, bj] * Xhat0[j]
            end
            if U != nothing
                arr[arr_length] = Whatk_blocks[i]
            end
            Xhatk[i] = (output_function == nothing) ?
                overapproximate(blocks[i], MinkowskiSumArray(arr)) :
                MinkowskiSumArray(copy(arr))
        end
        array = CartesianProductArray(copy(Xhatk))

        X_store = (output_function == nothing) ?
            array :
            box_approximation(output_function(array))

        terminate, skip, rs = termination(k, X_store, t0)
        if  !(isdisjoint(X_store, UnionSetArray(guards_proj)))
            for i in 1:bd
               bi = partition[diff_blocks[i]]
               for (j, bj) in enumerate(partition)
                   arr[j] = ϕpowerk[bi, bj] * Xhat0[j]
               end
               if U != nothing
                   arr[arr_length] = Whatk_diff_blocks[i]
               end
               Xhatk_d[i] = (output_function == nothing) ?
                   overapproximate(diff_blocks[i], MinkowskiSumArray(arr)) :
                   MinkowskiSumArray(copy(arr))
            end

            array_d = CartesianProductArray(copy(Xhatk_d))

            X_store_d = (output_function == nothing) ?
               array_d :
               box_approximation(output_function(array_d))
            if !terminate
               rs_oa = LazySets.Approximations.overapproximate(rs, CartesianProductArray, block_options)
               X_store = combine_cpas(rs_oa, X_store_d, blocks, diff_blocks)
            end
        end

        t0 = t1
        t1 += δ
        store!(res, k, X_store, t0, t1, N)

        if terminate
            break
        end

        if U != nothing
            for i in 1:b
                bi = partition[blocks[i]]
                Whatk_blocks[i] = overapproximate_inputs(k, blocks[i],
                    Whatk_blocks[i] + row(ϕpowerk, bi) * inputs)
            end
        end

         ##for diffs
        if U != nothing
            for i in 1:bd
                bi = partition[diff_blocks[i]]
                Whatk_diff_blocks[i] = overapproximate_inputs(k, diff_blocks[i],
                    Whatk_diff_blocks[i] + row(ϕpowerk, bi) * inputs)
            end
        end

         _A_mul_B!(ϕpowerk_cache, ϕpowerk, ϕ)
        copyto!(ϕpowerk, ϕpowerk_cache)

         k += 1
    end

    return k, skip
end
