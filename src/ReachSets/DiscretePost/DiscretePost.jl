"""
    DiscretePost

Abstract supertype of all discrete post operators.

### Notes

All discrete post operators should provide the following method, in addition
to those provided for general post operators:
```julia
tube⋂inv(op::DiscretePost, reach_tube, invariant, Rsets, start_interval)
```
"""
abstract type DiscretePost <: PostOperator end

function cluster(op::DiscretePost, reach_sets)
    # TODO apply some clustering
    return reach_sets
end

function isfixpoint(op::DiscretePost,
                    reach_set::ReachSet{LazySet{N}, N},
                    passed_list,
                    loc_id
                   ) where N
    if isassigned(passed_list, loc_id)
        for other_reach_set in passed_list[loc_id]
            if reach_set.X ⊆ other_reach_set.X
                info("found a fixpoint in some reach tube")
                return true
            end
        end
        return false
    else
        passed_list[loc_id] = Vector{ReachSet{LazySet{N}, N}}()
        return false
    end
end
