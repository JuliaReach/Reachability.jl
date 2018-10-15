import LazySets.use_precise_ρ

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

function cluster(op::DiscretePost,
                 reach_sets::Vector{ReachSet{LazySet{N}, N}},
                 options::Options) where N<:Real
    strategy = options[:clustering]
    if strategy == :none
        # no clustering
        return reach_sets
    elseif strategy == :chull
        # cluster all sets in a convex hull and overapproximate that set with
        # oct directions
        chull = ConvexHullArray(
            LazySet{N}[reach_set.X for reach_set in reach_sets])
        chull_oa = overapproximate(chull,
                                   Approximations.OctDirections(dim(chull)))
        return [ReachSet{LazySet{N}, N}(chull_oa, reach_sets[1].t_start,
                reach_sets[end].t_end)]
    end
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

# default: always apply line search
function use_precise_ρ(op::DiscretePost,
                             cap::Intersection{N})::Bool where N<:Real
    return true
end

function get_overapproximation_option(op::DiscretePost, n::Int)
    oa = op.options.dict[:overapproximation]
    if oa isa Symbol
        dirs = Utils.interpret_template_direction_symbol(oa)
        return dirs(n)
    elseif oa <: LazySets.LazySet
        return oa
    else
        error("received unknown :overapproximation option $oa")
    end
end
