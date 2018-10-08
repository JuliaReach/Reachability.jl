# ====================================================
# Textbook implementation of a discrete post operator.
# Uses concrete intersection.
#
### Notes
#
# The current implementation requires that the `Polyhedra` # library is loaded,
# because some concrete operations between polytopes are used.
#
# Currently, we assume that source invariants, target invariants and guards are
# polytopes in constraint representation resp. HalfSpaces.
#
### Algorithm
#
# The algorithm is based on [Flowpipe-Guard Intersection for Reachability
# Computations with Support Functions](
# http://spaceex.imag.fr/sites/default/files/frehser_adhs2012.pdf).
# ====================================================

struct TextbookDiscretePost <: DiscretePost
end

function init(op::TextbookDiscretePost, system, options_input)
    @assert isdefined(Main, :Polyhedra) "this algorithm needs the package " *
            "'Polyhedra' to be loaded"

    options_input.dict[:n] = statedim(system, 1)

    # solver-specific options (adds default values for unspecified options)
    options = validate_solver_options_and_add_default_values!(options_input)

    # Input -> Output variable mapping
    options.dict[:inout_map] =
        inout_map_reach(options[:partition], options[:blocks], options[:n])

    return options
end

function tube⋂inv!(op::TextbookDiscretePost,
                   reach_tube::Vector{<:ReachSet{<:LazySet{N}}},
                   invariant,
                   Rsets,
                   start_interval
                  ) where {N}
    # take intersection with source invariant

    if invariant isa HalfSpace
        # TODO temporary conversion to HPolytope
        invariant = HPolytope([invariant])
    else
        @assert invariant isa HPolytope
    end

    # TODO First check for empty intersection, which can be more efficient.
    #      However, we need to make sure that the emptiness check does not just
    #      compute the concrete intersection; otherwise, we would do the work
    #      twice. This is currently the case for 'Polyhedra' polytopes.
    intersections = Vector{ReachSet{LazySet{N}, N}}()
    for reach_set in reach_tube
        rs = reach_set.X
        # TODO temporary workaround for 1D sets
        if dim(rs) == 1
            reach_tube = VPolytope(vertices_list(
                Approximations.overapproximate(rs, LazySets.Interval)))
        # TODO offer a lazy intersection here
        # TODO offer more options instead of taking the VPolytope intersection
        elseif rs isa CartesianProductArray
            reach_tube = VPolytope(vertices_list(rs))
        else
            error("unsupported set type for reach tube: $(typeof(rs))")
        end
        intersect = intersection(invariant, reach_tube)
        if isempty(intersect)
            break
        end
        push!(intersections, ReachSet{LazySet{N}, N}(intersect,
            reach_set.t_start + start_interval[1],
            reach_set.t_end + start_interval[2]))
    end

    append!(Rsets, intersections)
end

function post(op::TextbookDiscretePost,
              system::InitialValueProblem{<:HybridSystem, <:LazySet{N}},
              waiting_list,
              passed_list,
              source_loc_id,
              tube⋂inv,
              jumps
             ) where {N}
    HS = system.s
    jumps += 1
    for trans in out_transitions(HS, source_loc_id)
        info("Considering transition: $trans")
        target_loc_id = target(HS, trans)
        target_loc = HS.modes[target(HS, trans)]
        target_invariant = target_loc.X
        trans_annot = HS.resetmaps[symbol(HS, trans)]
        guard = trans_annot.X
        assignment = trans_annot.A

        # TODO temporary conversion to HPolytope
        @assert target_invariant isa HalfSpace
        target_invariant = HPolytope([target_invariant])
        @assert guard isa HPolytope

        # perform jumps
        post_jump = Vector{ReachSet{LazySet{N}, N}}()
        sizehint!(post_jump, length(tube⋂inv))
        for reach_set in tube⋂inv
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
        clustered = cluster(op, post_jump)

        # push new sets after jump (unless a fixpoint is detected)
        for reach_set in clustered
            if passed_list != nothing
                if isfixpoint(op, reach_set, passed_list, target_loc_id)
                    continue
                end
                push!(passed_list[target_loc_id], reach_set)
            end
            push!(waiting_list, (target(HS, trans), reach_set, jumps))
        end
    end
end

function cluster(op::TextbookDiscretePost, reach_sets)
    # TODO apply some clustering
    return reach_sets
end

function isfixpoint(op::TextbookDiscretePost,
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
