function compute_dimensions(partition, blocks)
    dimensions = Vector{Int}()
    dims = 0
    next_idx = 1
    for block_idx in blocks
        while next_idx < block_idx
            # sum up dimensions of skipped blocks
            dims += length(partition[next_idx])
            next_idx += 1
        end
        append!(dimensions, partition[next_idx])
        next_idx += 1
    end
    return dimensions
end
