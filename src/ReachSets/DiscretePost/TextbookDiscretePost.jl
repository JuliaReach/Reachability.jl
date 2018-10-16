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
    options::Options
end

function TextbookDiscretePost()
    defaults = Options()
    setindex!(defaults, Hyperrectangle, :overapproximation)
    return TextbookDiscretePost(defaults)
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
                  )::Vector{ReachSet{LazySet{N}, N}} where {N}
    # take intersection with source invariant

    # TODO First check for empty intersection, which can be more efficient.
    #      However, we need to make sure that the emptiness check does not just
    #      compute the concrete intersection; otherwise, we would do the work
    #      twice. This is currently the case for 'Polyhedra' polytopes.
    intersections = Vector{ReachSet{LazySet{N}, N}}()
    for reach_set in reach_tube
        rs = reach_set.X
        @assert rs isa CartesianProductArray
        if length(array(rs)) == 1
            # TODO workaround for lazy X0
            rs_converted = Approximations.overapproximate(rs,
                Approximations.BoxDirections(dim(rs)))
        else
            rs_converted = HPolytope(constraints_list(rs))
        end
        R⋂I = intersection(invariant, rs_converted)
        if isempty(R⋂I)
            break
        end
        push!(intersections, ReachSet{LazySet{N}, N}(R⋂I,
            reach_set.t_start + start_interval[1],
            reach_set.t_end + start_interval[2]))
    end

    append!(Rsets, intersections)
    return intersections
end

function post(op::TextbookDiscretePost,
              HS::HybridSystem,
              waiting_list::Vector{Tuple{Int, ReachSet{LazySet{N}, N}, Int}},
              passed_list,
              source_loc_id,
              tube⋂inv,
              jumps,
              options
             ) where {N}
    jumps += 1
    for trans in out_transitions(HS, source_loc_id)
        info("Considering transition: $trans")
        target_loc_id = target(HS, trans)
        target_loc = HS.modes[target(HS, trans)]
        target_invariant = target_loc.X
        trans_annot = HS.resetmaps[symbol(HS, trans)]
        guard = trans_annot.X
        assignment = trans_annot.A

        # perform jumps
        post_jump = Vector{ReachSet{LazySet{N}, N}}()
        sizehint!(post_jump, length(tube⋂inv))
        for reach_set in tube⋂inv
            # check intersection with guard
            R⋂G = intersection(guard, HPolytope(constraints_list(reach_set.X)))
            if isempty(R⋂G)
                continue
            end
            # apply assignment
            # TODO converting to HPolytope; this should be handled automatically
            A⌜R⋂G⌟ = linear_map(assignment, R⋂G, output_type=HPolytope)
            # intersect with target invariant
            A⌜R⋂G⌟⋂I = intersection(target_invariant, A⌜R⋂G⌟)
            if isempty(A⌜R⋂G⌟⋂I)
                continue
            end

            # store result
            push!(post_jump, ReachSet{LazySet{N}, N}(A⌜R⋂G⌟⋂I,
                                                     reach_set.t_start,
                                                     reach_set.t_end))
        end

        postprocess(op, HS, post_jump, options, waiting_list, passed_list,
            target_loc_id, jumps)
    end
end
