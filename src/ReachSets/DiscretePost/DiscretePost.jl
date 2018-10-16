import LazySets.use_precise_ρ

"""
    DiscretePost

Abstract supertype of all discrete post operators.

### Notes

All discrete post operators should provide the following method, in addition
to those provided for general post operators:
```julia
tube⋂inv!(op::DiscretePost, reach_tube::Vector{<:ReachSet{<:LazySet{N}}},
          invariant, Rsets, start_interval)::Vector{ReachSet{LazySet{N}, N}}
```
"""
abstract type DiscretePost <: PostOperator end

function postprocess(op,
                     HS,
                     post_jump,
                     options,
                     waiting_list,
                     passed_list,
                     target_loc_id,
                     jumps)
    fixpoint_strategy = options[:fixpoint_check]

    if fixpoint_strategy == :eager
        # eager fixpoint checking
        post_jump_filtered =
            filter(x -> !isfixpoint(op, x, passed_list, target_loc_id),
                   post_jump)
    else
        post_jump_filtered = post_jump
    end

    if (isempty(post_jump_filtered))
        # fixpoint found or transition can never be taken
        return
    end

    # apply clustering
    clustered = cluster(op, post_jump_filtered, options)

    # push new sets after jump (unless a fixpoint is detected)
    for reach_set in clustered
        if fixpoint_strategy != :none
            if fixpoint_strategy == :lazy &&
                    isfixpoint(op, reach_set, passed_list, target_loc_id)
                continue
            end
            push!(passed_list[target_loc_id], reach_set)
        end
        push!(waiting_list, (target_loc_id, reach_set, jumps))
    end
end

function cluster(op::DiscretePost,
                 reach_sets::Vector{ReachSet{LazySet{N}, N}},
                 options::Options) where N<:Real
    clustering_strategy = options[:clustering]
    dirs = op.options.[:overapproximation]
    if clustering_strategy == :none
        # no clustering
        return reach_sets
    elseif clustering_strategy == :chull
        # cluster all sets in a convex hull and overapproximate that set with
        # oct directions
        chull = ConvexHullArray(
            LazySet{N}[reach_set.X for reach_set in reach_sets])
        chull_oa = overapproximate(chull,
                                   dirs)
        return [ReachSet{LazySet{N}, N}(chull_oa, reach_sets[1].t_start,
                reach_sets[end].t_end)]
    end
end

function isfixpoint(op::DiscretePost,
                    reach_set::ReachSet{LazySet{N}, N},
                    passed_list,
                    loc_id
                   ) where N
    @assert passed_list != nothing
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
