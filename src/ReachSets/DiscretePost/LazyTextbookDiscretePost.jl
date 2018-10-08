# ==============================================================================
# Textbook implementation of a discrete post operator, but with lazy operations.
# ==============================================================================

struct LazyTextbookDiscretePost <: DiscretePost
    options::Options
end

LazyTextbookDiscretePost() = LazyTextbookDiscretePost(Options())

function init(op::LazyTextbookDiscretePost, system, options_input)
    options_input.dict[:n] = statedim(system, 1)

    # solver-specific options (adds default values for unspecified options)
    options = validate_solver_options_and_add_default_values!(options_input)

    # Input -> Output variable mapping
    options.dict[:inout_map] =
        inout_map_reach(options[:partition], options[:blocks], options[:n])

    # set up operator-specific options
    # TODO temporarily use box approximation
    @assert !haskey(op.options.dict, :overapproximation)
    op.options.dict[:overapproximation] = Hyperrectangle

    return options
end

function tube⋂inv!(op::LazyTextbookDiscretePost,
                   reach_tube::Vector{<:ReachSet{<:LazySet{N}}},
                   invariant,
                   Rsets,
                   start_interval
                  ) where {N}
    intersections = Vector{ReachSet{LazySet{N}, N}}()
    for reach_set in reach_tube
        R⋂I = Intersection(invariant, reach_set.X)
        if isempty(R⋂I)
            break
        end
        push!(intersections, ReachSet{LazySet{N}, N}(R⋂I,
            reach_set.t_start + start_interval[1],
            reach_set.t_end + start_interval[2]))
    end

    append!(Rsets, intersections)
end

function post(op::LazyTextbookDiscretePost,
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

        # perform jumps
        post_jump = Vector{ReachSet{LazySet{N}, N}}()
        sizehint!(post_jump, length(tube⋂inv))
        for reach_set in tube⋂inv
            # check intersection with guard
            R⋂G = Intersection(reach_set.X, guard)
            if isempty(R⋂G)
                continue
            end
            # apply assignment
            A⌜R⋂G⌟ = LinearMap(assignment, R⋂G)
            # intersect with target invariant
            A⌜R⋂G⌟⋂I = Intersection(target_invariant, A⌜R⋂G⌟)
            # overapproximate final set
            res = overapproximate(A⌜R⋂G⌟⋂I, op.options[:overapproximation];
                                  upper_bound=true)

            # check if the final set is empty
            if isempty(res)
                continue
            end

            # store result
            push!(post_jump, ReachSet{LazySet{N}, N}(res,
                                                     reach_set.t_start,
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
            if !isfixpoint(op, reach_set, passed_list, target_loc_id)
                push!(passed_list[target_loc_id], reach_set)
                push!(waiting_list, (target(HS, trans), reach_set, jumps))
            end
        end
    end
end
