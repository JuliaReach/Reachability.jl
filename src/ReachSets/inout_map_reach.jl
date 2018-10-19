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
    inout = Vector{Int}(undef, n)
    blocks_idx = 1
    next_block_idx = blocks[blocks_idx]
    out_idx = 1
    for (i, bi) in enumerate(partition)
        if next_block_idx == i
            # block/variables present in output
            inout[bi] = out_idx:(out_idx + length(bi) - 1)
            out_idx += length(bi)
            if blocks_idx < length(blocks)
                blocks_idx += 1
                next_block_idx = blocks[blocks_idx]
            else
                next_block_idx = 0
            end
        else
            @assert next_block_idx > i || next_block_idx == 0
            # block/variables not present in output
            inout[bi] .= 0
        end
    end
    return inout
end
