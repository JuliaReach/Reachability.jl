import LazySets.use_precise_ρ

"""
    DiscretePost

Abstract supertype of all discrete post operators.

### Notes

All discrete post operators should provide the following method, in addition
to those provided for general post operators:
```julia
tube⋂inv!(𝒫::DiscretePost, reach_tube::Vector{<:ReachSet{<:LazySet{N}}},
          invariant, Rsets, start_interval)::Vector{ReachSet{LazySet{N}, N}}
```
"""
abstract type DiscretePost <: PostOperator end

function postprocess(𝒫,
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
            filter(x -> !isfixpoint(𝒫, x, passed_list, target_loc_id),
                   post_jump)
    else
        post_jump_filtered = post_jump
    end

    if (isempty(post_jump_filtered))
        # fixpoint found or transition can never be taken
        return
    end

    # apply clustering
    clustered = cluster(𝒫, post_jump_filtered, options)

    # push new sets after jump (unless a fixpoint is detected)
    for reach_set in clustered
        if fixpoint_strategy != :none
            if fixpoint_strategy == :lazy &&
                    isfixpoint(𝒫, reach_set, passed_list, target_loc_id)
                continue
            end
            push!(passed_list[target_loc_id], reach_set)
        end
        push!(waiting_list, (target_loc_id, reach_set, jumps))
    end
end

function cluster(𝒫::DiscretePost,
                 reach_sets::Vector{ReachSet{LazySet{N}, N}},
                 options::Options) where N<:Real
    clustering_strategy = options[:clustering]
    dirs = 𝒫.options[:overapproximation]
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

function isfixpoint(𝒫::DiscretePost,
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
function use_precise_ρ(𝒫::DiscretePost,
                             cap::Intersection{N})::Bool where N<:Real
    return true
end

function get_overapproximation_option(𝒫::DiscretePost, n::Int)
    oa = 𝒫.options[:overapproximation]
    if oa isa Symbol
        dirs = Utils.interpret_template_direction_symbol(oa)
        return dirs(n)
    elseif oa <: LazySets.LazySet
        return oa
    else
        error("received unknown :overapproximation option $oa")
    end
end

# --- default methods for handling assignments ---

# default implementation: use 'apply' from MathematicalSystems
function apply_assignment(𝒫::DiscretePost,
                          constrained_map::AbstractMap,
                          R⋂G::LazySet;
                          kwargs...)
    return apply(constrained_map, R⋂G)
end

# for reset maps: return a lazy ResetMap from LazySets
function apply_assignment(𝒫::DiscretePost,
                          constrained_map::ConstrainedResetMap,
                          R⋂G::LazySet;
                          kwargs...)
    return LazySets.ResetMap(R⋂G, constrained_map.dict)
end
