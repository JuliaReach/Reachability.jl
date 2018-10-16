# ==============================================================================
# Textbook implementation of a discrete post operator, but with lazy operations.
# ==============================================================================

import LazySets.use_precise_ρ

struct LazyTextbookDiscretePost <: DiscretePost
    options::Options
end

# default options for the LazyTextbookDiscretePost discrete post operator
function LazyTextbookDiscretePost()
    defaults = Options()
    setindex!(defaults, Hyperrectangle, :overapproximation)
    setindex!(defaults, false, :lazy_R⋂I)
    setindex!(defaults, true, :lazy_R⋂G)
    setindex!(defaults, true, :lazy_A⌜R⋂G⌟⋂I)
    return LazyTextbookDiscretePost(defaults)
end

function init(op::LazyTextbookDiscretePost, system, options_input)
    options_input.dict[:n] = statedim(system, 1)

    # solver-specific options (adds default values for unspecified options)
    options = validate_solver_options_and_add_default_values!(options_input)

    # Input -> Output variable mapping
    options.dict[:inout_map] =
        inout_map_reach(options[:partition], options[:blocks], options[:n])

    # set up operator-specific options
    @assert haskey(op.options.dict, :overapproximation)
    

    return options
end

function tube⋂inv!(op::LazyTextbookDiscretePost,
                   reach_tube::Vector{<:ReachSet{<:LazySet{N}}},
                   invariant,
                   Rsets,
                   start_interval
                  )::Vector{ReachSet{LazySet{N}, N}} where {N}
    intersections = Vector{ReachSet{LazySet{N}, N}}()
    dirs = op.options[:overapproximation]
    for reach_set in reach_tube
        R⋂I = Intersection(invariant, reach_set.X)
        if isempty(R⋂I)
            break
        end
        if !op.options[:lazy_R⋂I]
            R⋂I = overapproximate(R⋂I, dirs)
        end
        push!(intersections, ReachSet{LazySet{N}, N}(R⋂I,
            reach_set.t_start + start_interval[1],
            reach_set.t_end + start_interval[2]))
    end

    append!(Rsets, intersections)
    return intersections
end

function post(op::LazyTextbookDiscretePost,
              HS::HybridSystem,
              waiting_list::Vector{Tuple{Int, ReachSet{LazySet{N}, N}, Int}},
              passed_list,
              source_loc_id,
              tube⋂inv,
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
        sizehint!(post_jump, length(tube⋂inv))
        for reach_set in tube⋂inv
            # check intersection with guard
            R⋂G = Intersection(reach_set.X, guard)
            if isempty(R⋂G)
                continue
            end
            if !op.options[:lazy_R⋂G]
               R⋂G = overapproximate(R⋂G, dirs)
            end
            # apply assignment
            A⌜R⋂G⌟ = LinearMap(assignment, R⋂G)
            # intersect with target invariant
            A⌜R⋂G⌟⋂I = Intersection(target_invariant, A⌜R⋂G⌟)
            # overapproximate final set
            if !op.options[:lazy_A⌜R⋂G⌟⋂I]
                res = overapproximate(A⌜R⋂G⌟⋂I, dirs)
            else
                res = A⌜R⋂G⌟⋂I
            end

            # check if the final set is empty
            if isempty(res)
                continue
            end

            # store result
            push!(post_jump, ReachSet{LazySet{N}, N}(res,
                                                     reach_set.t_start,
                                                     reach_set.t_end))
        end

        postprocess(op, HS, post_jump, options, waiting_list, passed_list,
            target_loc_id, jumps)
    end
end

# --- line search policies ---

# usually do not use line search
function use_precise_ρ(op::LazyTextbookDiscretePost,
                             cap::Intersection{N})::Bool where N<:Real
    return false
end

# use line search for the outermost level, which is a LinearMap
function use_precise_ρ(op::LazyTextbookDiscretePost,
                             cap::Intersection{N, <:LinearMap{N}}
                            )::Bool where N<:Real
    return true
end
