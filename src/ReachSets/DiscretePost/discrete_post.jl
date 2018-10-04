function discrete_post!(waiting_list, passed_list, HS, cur_loc_id,
                        reach_tube_in_invariant, jumps, N)
    jumps += 1
    for trans in out_transitions(HS, cur_loc_id)
        info("Considering transition: $trans")
        target_loc_id = target(HS, trans)
        target_loc = HS.modes[target(HS, trans)]
        target_invariant = target_loc.X
        guard = HS.resetmaps[target_loc_id].X
        assignment = HS.resetmaps[target_loc_id].A

        # TODO temporary conversion to HPolytope
        @assert target_invariant isa HalfSpace
        target_invariant = HPolytope([target_invariant])
        @assert guard isa HPolytope

        # perform jumps
        post_jump = Vector{ReachSet{LazySet{N}, N}}()
        sizehint!(post_jump, length(reach_tube_in_invariant))
        for reach_set in reach_tube_in_invariant
            # check intersection with guard
            cap = intersection(guard, VPolytope(vertices_list(reach_set.X)))
            if isempty(cap)
                continue
            end
            # apply assignment
            asgn = linear_map(assignment, cap)
            # intersect with target invariant
            res = intersection(target_invariant, asgn)
            if isempty(res)
                continue
            end

            # store result
            push!(post_jump, ReachSet{LazySet{N}, N}(res, reach_set.t_start,
                                                      reach_set.t_end))
        end

        if (isempty(post_jump))
            # transition can never be taken
            continue
        end

        # apply clustering
        clustered = cluster(post_jump)

        # push new sets after jump (unless a fixpoint is detected)
        for reach_set in clustered
            if !isfixpoint(reach_set, passed_list, target_loc_id)
                push!(passed_list[target_loc_id], reach_set)
                push!(waiting_list, (target(HS, trans), reach_set, jumps))
            end
        end
    end
end

function cluster(sets)
    # TODO apply some clustering
    return sets
end

function isfixpoint(reach_set::ReachSet{LazySet{N}, N},
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
