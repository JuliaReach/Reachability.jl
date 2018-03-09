"""
    inout_map_reach(partition::AbstractVector{<:AbstractVector{Int}},
                    blocks::AbstractVector{Int},
                    n::Int)

Compute a mapping from input dimension to output dimension.
The map returns `0` for all input dimensions not present in the output.

### Input

- `partition` -- block partition; elements are start and end indices of a block
- `blocks`    -- list of all output block indices in the partition
- `n`         -- total number of input dimensions

### Output

A vector mapping from input dimension (index) to output dimension (entry).
"""
function inout_map_reach(partition::AbstractVector{<:AbstractVector{Int}},
                         blocks::AbstractVector{Int},
                         n::Int)
    inout = Vector{Int}(n)
    blocks_idx = 1
    out_idx = 1
    for (i, bi) in enumerate(partition)
        if blocks[blocks_idx] != i
            @assert blocks[blocks_idx] > i
            # block/variables not present in output
            inout[bi] = 0
        else
            # block/variables present in output
            inout[bi] = out_idx:(out_idx + length(bi) - 1)
            out_idx += length(bi)
            blocks_idx += 1
        end
    end
    return inout
end
