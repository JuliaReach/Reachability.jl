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

function tube⋂inv!(op::ApproximatingDiscretePost,
                   reach_tube::Vector{<:ReachSet{<:LazySet{N}}},
                   invariant,
                   Rsets,
                   start_interval
                  ) where {N}

    dirs = get_overapproximation_option(op, dim(invariant))

    # counts the number of sets R⋂I added to Rsets
    count = 0
    for reach_set in reach_tube
        R⋂I = Intersection(invariant, reach_set.X)
        if op.options[:check_invariant_intersection] && isempty(R⋂I)
            break
        end
        # return an overapproximation
        push!(Rsets, ReachSet{LazySet{N}, N}(
            overapproximate(R⋂I, dirs),
            reach_set.t_start + start_interval[1],
            reach_set.t_end + start_interval[2]))
        count = count + 1
    end
    return count
end

function post(op::ApproximatingDiscretePost,
              HS::HybridSystem,
              waiting_list::Vector{Tuple{Int, ReachSet{LazySet{N}, N}, Int}},
              passed_list,
              source_loc_id,
              tube⋂inv,
              count_Rsets,
              jumps,
              options
             ) where {N}
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
        sizehint!(post_jump, count_Rsets)
        for reach_set in tube⋂inv[length(tube⋂inv) - count_Rsets + 1 : end]
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
end