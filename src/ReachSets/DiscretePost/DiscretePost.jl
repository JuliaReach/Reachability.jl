"""
    DiscretePost

Abstract supertype of all discrete post operators.

### Notes

All discrete post operators should provide the following methods, in addition
to those provided for general post operators:
```julia
tube⋂inv(op::DiscretePost, reach_tube, invariant, Rsets, start_interval)
cluster(op::DiscretePost, reach_sets)
isfixpoint(op::DiscretePost, reach_set, passed_list, loc_id)
```
"""
abstract type DiscretePost <: PostOperator end
