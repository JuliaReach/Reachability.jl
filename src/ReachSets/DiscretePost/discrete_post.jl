function discrete_post!(waiting_list, HS, cur_loc_id, reach_tube_in_invariant,
                        jumps, N)
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

        # push new sets after jump
        for set in clustered
            push!(waiting_list, (target(HS, trans), set, jumps + 1))
        end
    end
end

function cluster(sets)
    # TODO apply some clustering
    return sets
end
