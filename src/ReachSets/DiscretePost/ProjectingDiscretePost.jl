# ==============================================================================
# Implementation of a discrete post operator that projects guards to constrained
# dimensions and computes the reach tubes only in low dimensions as well.
# ==============================================================================

mutable struct ReachTubeAndInvariant
    R::Vector{<:ReachSet{<:LazySet{N}}}
    I::LazySet
    Rsets::Vector{ReachSet{LazySet{N}, N}}
    start_interval::Vector{Int}
end

struct ProjectingDiscretePost <: DiscretePost
    options::Options
    RI::Union{ReachTubeAndInvariant, Nothing}
end

ProjectingDiscretePost() = ProjectingDiscretePost(Options(), nothing)

function init(op::ProjectingDiscretePost, system, options_input)
    options_input.dict[:n] = statedim(system, 1)

    # solver-specific options (adds default values for unspecified options)
    options = validate_solver_options_and_add_default_values!(options_input)

    # Input -> Output variable mapping
    options.dict[:inout_map] =
        inout_map_reach(options[:partition], options[:blocks], options[:n])

    # set up operator-specific options
    @assert !haskey(op.options.dict, :loc2vars)
    op.options.dict[:loc2vars] = find_constrained_dimensions(system)
    # TODO use hyperrectangles for now
    @assert !haskey(op.options.dict, :overapproximation)
    op.options.dict[:overapproximation] = Hyperrectangle

    return options
end

function tube⋂inv!(op::ProjectingDiscretePost,
                   reach_tube::Vector{<:ReachSet{<:LazySet{N}}},
                   invariant,
                   Rsets,
                   start_interval
                  ) where {N}
    # no-op, just store the data for later
    op.RI = ReachTubeAndInvariant(reach_tube, invariant, Rsets, start_interval)
end

function tube⋂inv_delayed!(op::ProjectingDiscretePost,
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

function post(op::ProjectingDiscretePost,
              HS::HybridSystem,
              waiting_list::Vector{Tuple{Int, ReachSet{LazySet{N}, N}, Int}},
              passed_list,
              source_loc_id,
              tube⋂inv,
              jumps
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

function find_constrained_dimensions(HS)
    n = statedim(HS, 1)
    m = length(HS.modes)
    loc2vars = Vector{AbstractVector{Int}}(undef, m)
    for loc_id in states(HS)
        relevant_bitmap = Vector{Bool}(undef, n)
        all_relevant = false
        has_transitions = false
        for trans in out_transitions(HS, loc_id)
            has_transitions = true
            guard = HS.resetmaps[symbol(HS, trans)].X
            if guard isa Union{HalfSpace, Hyperplane, Line}
                for i in constraint_dimensions(guard)
                    relevant_bitmap[i] = true
                end
            else
                all_relevant = true
                break
            end
        end
        if !has_transitions
            # just choose the first variable
            relevant = [1]
        elseif all_relevant
            relevant = 1:n
            info("considering *all* dimensions in location $loc_id")
        else
            relevant = Int[]
            for i in 1:length(relevant_bitmap)
                if relevant_bitmap[i]
                    push!(relevant, i)
                end
            end
            @assert !isempty(relevant) "did not find any relevant dimensions " *
                "in location $loc_id"
            info("considering dimensions $relevant in location $loc_id")
        end
        loc2vars[loc_id] = relevant
    end
    return loc2vars
end
