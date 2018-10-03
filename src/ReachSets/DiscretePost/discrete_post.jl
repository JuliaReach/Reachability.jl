function discrete_post!(waiting_list, HS, cur_loc_id, intersectedRset, jumps)
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

        # check intersection with guard
        rsetIntersMinus = [intersection(guard, VPolytope(vertices_list(hi))) for hi in intersectedRset]
        filter!(!isempty, rsetIntersMinus)
        if (isempty(rsetIntersMinus))
            continue
        end

        # apply assignment
        rsetIntersMinus = [linear_map(assignment, ri) for ri in rsetIntersMinus]

        # check intersection with target invariant
        info("Intersection with I\^+")
        rsetIntersPlus = [intersection(target_invariant, hi) for hi in rsetIntersMinus]
        filter!(!isempty, rsetIntersPlus)

        # push new sets after jump
        rsetHull = ConvexHullArray(rsetIntersPlus)
        for rh in rsetIntersPlus
            push!(waiting_list, (target(HS, trans), rh, jumps + 1))
        end

        # TODO check intersection with forbidden states
        # push!(waiting_list,(target(HS, trans), rsetHull))
    end
end
