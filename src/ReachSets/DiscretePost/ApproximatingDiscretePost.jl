# ==================================================================
# Textbook implementation of a discrete post operator, but with lazy
# intersection followed by an overapproximation.
# ==================================================================

struct ApproximatingDiscretePost <: DiscretePost
    options::Options
end

function ApproximatingDiscretePost()
    defaults = Options()
    setindex!(defaults, Hyperrectangle, :overapproximation)
    setindex!(defaults, false, :check_invariant_intersection)
    return ApproximatingDiscretePost(defaults)
end

function init(op::ApproximatingDiscretePost, system, options_input)
    options_input.dict[:n] = statedim(system, 1)

    # solver-specific options (adds default values for unspecified options)
    options = validate_solver_options_and_add_default_values!(options_input)

    # Input -> Output variable mapping
    options.dict[:inout_map] =
        inout_map_reach(options[:partition], options[:blocks], options[:n])

    # check operator-specific options
    @assert haskey(op.options.dict, :overapproximation)

    return options
end

function post(op::ApproximatingDiscretePost,
              HS::HybridSystem,
              waiting_list::Vector{Tuple{Int, ReachSet{LazySet{N}, N}, Int}},
              passed_list,
              reach_tube::Vector{<:ReachSet{<:LazySet{N}}},
              Rsets,
              start_interval,
              source_loc_id,
              jumps,
              options
             )::Void where {N}
    # intersection with source invariant
    source_invariant = HS.modes[source_loc_id].X
    tube⋂inv = Vector{ReachSet{LazySet{N}, N}}()
    dirs = get_overapproximation_option(op, dim(source_invariant))
    for reach_set in reach_tube
        R⋂I = Intersection(reach_set.X, source_invariant)
        if op.options[:check_invariant_intersection] && isempty(R⋂I)
            break
        end
        # return an overapproximation
        push!(tube⋂inv, ReachSet{LazySet{N}, N}(
            overapproximate(R⋂I, dirs),
            reach_set.t_start + start_interval[1],
            reach_set.t_end + start_interval[2]))
    end

    append!(Rsets, tube⋂inv)

    # terminate if no jump should be taken anymore
    if jumps == options[:max_jumps]
        return nothing
    end

    jumps += 1
    dirs = get_overapproximation_option(op, options[:n])
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
            A⌜R⋂G⌟o = overapproximate(A⌜R⋂G⌟, dirs)

            # intersect with target invariant
            A⌜R⋂G⌟o⋂I = Intersection(target_invariant, A⌜R⋂G⌟o)

            # check if the final set is empty
            if isempty(A⌜R⋂G⌟o⋂I)
                continue
            end

            # overapproximate final set once more
            res = overapproximate(A⌜R⋂G⌟o⋂I, dirs)

            # store result
            push!(post_jump, ReachSet{LazySet{N}, N}(res,
                                                     reach_set.t_start,
                                                     reach_set.t_end))
        end

        postprocess(op, HS, post_jump, options, waiting_list, passed_list,
            target_loc_id, jumps)
    end
    return nothing
end
