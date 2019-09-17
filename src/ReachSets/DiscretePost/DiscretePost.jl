import LazySets.use_precise_ρ

"""
    DiscretePost

Abstract supertype of all discrete post operators.

### Notes

All discrete post operators should provide the following method, in addition
to those provided for general post operators:
```julia
tube⋂inv!(𝒫::DiscretePost, reach_tube::Vector{<:AbstractReachSet{<:LazySet, N}},
          invariant, Rsets, start_interval) where {N}
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
                 reach_sets::Vector{RSN},
                 options::Options) where {SN, RSN<:AbstractReachSet{SN}}
    clustering_strategy = options[:clustering]
    oa = 𝒫.options[:overapproximation]
    if clustering_strategy == :none
        # no clustering, keeping original set
        return reach_sets
    elseif clustering_strategy == :none_oa
        # no clustering but overapproximation
        return [RSN(overapproximate(set(reach_set), oa),
        time_start(reach_set), time_end(reach_set)) for reach_set in reach_sets]
    elseif clustering_strategy == :chull
        # cluster all sets in a convex hull and overapproximate that set
        chull = ConvexHullArray([set(reach_set) for reach_set in reach_sets])
        chull_oa = overapproximate(chull, oa)
        return [RSN(chull_oa, time_start(reach_sets[1]),
                    time_end(reach_sets[end]))]
    end
end

function isfixpoint(𝒫::DiscretePost,
                    reach_set::RSN,
                    passed_list,
                    loc_id
                   ) where {SN, RSN<:AbstractReachSet{SN}}
    @assert passed_list != nothing
    if isassigned(passed_list, loc_id)
        for other_reach_set in passed_list[loc_id]
            if set(reach_set) ⊆ set(other_reach_set)
                info("found a fixpoint in some reach tube")
                return true
            end
        end
        return false
    else
        passed_list[loc_id] = Vector{RSN}()
        return false
    end
end

# default: always apply line search
function use_precise_ρ(𝒫::DiscretePost,
                             cap::Intersection{N})::Bool where N<:Real
    return true
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
