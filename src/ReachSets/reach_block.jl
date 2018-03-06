#=
    @reach_block!(ϕ, Xhat0, U, n, b, N, overapproximate, blocks, res)

Reachability computation of a given number of two-dimensional blocks of an
affine system with nondeterministic inputs.

The variants have the following structure:

### Input

- `ϕ`      -- sparse matrix of a discrete affine system
- `Xhat0`  -- initial set as a cartesian product over 2d blocks
- `U`      -- input set of undeterministic inputs
- `n`      -- ambient dimension
- `b`      -- number of blocks
- `N`      -- number of sets computed
- `overapproximate` -- function for overapproximation
- `blocks` -- the block indices to be computed
- `res`    -- storage space for the result, a linear array of CartesianProductArray

### Output

Array of the cartesian product of two-dimensional sets for the given block
indices, and ZeroSet's for the rest of them.
It is obtained by reachability computation of a discrete affine system with
nondeterministic inputs.
=#

# sparse, with input
macro reach_block_sparse!(SINGLE, BLOCKS_ARG, STORE)
     @eval begin
         function reach_block!(ϕ::SparseMatrixCSC{NUM, Int},
                                       Xhat0::Vector{<:LazySet{NUM}},
                                       U::ConstantNonDeterministicInput,
                                       n::Int,
                                       b::Int,
                                       N::Int,
                                       overapproximate::Function,
                                       $BLOCKS_ARG,
                                       res::Vector{<:LazySet{NUM}}
                                      )::Void where {NUM}
            res[1] = $SINGLE ? Xhat0[bi] :
                               CartesianProductArray(Xhat0)
            if N == 1
                return nothing
            end

            @inline F(bi::Int, bj::Int) = ϕpowerk[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]
            @inline oa_inputs_1(bi::Int) =
                overapproximate(sparse(1:2, (2*bi-1):(2*bi), [1., 1.], 2, n) * inputs)
            @inline oa_inputs_2(bi::Int, Whatk_bi) =
                overapproximate(Whatk_bi + ϕpowerk[(2*bi-1):(2*bi), :] * inputs)

            ϕpowerk = copy(ϕ)
            dummy_set = ZeroSet(2)
            inputs = next_set(U)

            if $SINGLE
                Whatk_bi = oa_inputs_1(bi)
            else
                Xhatk = Vector{LazySet{NUM}}(b)
                @inbounds for bi in 1:b
                     Xhatk[bi] = dummy_set
                end
                Whatk = Vector{LazySet{NUM}}(b)
                @inbounds for bi in blocks
                    Whatk[bi] = oa_inputs_1(bi)
                end
            end

            k = 2
            @inbounds while true
                @inbounds for bi in blocks
                    Xhatk_bi = dummy_set
                    @inbounds for bj in 1:b
                        if findfirst(F(bi, bj)) != 0
                            Xhatk_bi = Xhatk_bi + F(bi, bj) * Xhat0[bj]
                        end
                    end
                    $STORE = overapproximate(Xhatk_bi + Whatk[bi])
                end
                res[k] = $SINGLE ? $STORE :
                                   res[k] = CartesianProductArray(copy(Xhatk))

                if k == N
                    break
                end

                if $SINGLE
                    Whatk_bi = oa_inputs_2(bi, Whatk_bi)
                else
                    @inbounds for bi in blocks
                        Whatk[bi] = oa_inputs_2(bi, Whatk[bi])
                    end
                end

                ϕpowerk = ϕpowerk * ϕ
                k += 1
            end
        end
    end
    return nothing
end

# single block
@reach_block_sparse!(true, # single block?
                     bi::Int, # blocks argument
                     Xhatk[bi] # temporary storage
                    )
# multiple blocks
@reach_block_sparse!(false, # single block?
                     blocks::AbstractVector{Int}, # blocks argument
                     Xhatk_bi # temporary storage
                    )
