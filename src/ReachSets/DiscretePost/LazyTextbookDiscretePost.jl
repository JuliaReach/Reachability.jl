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
    setindex!(defaults, false, :check_invariant_intersection)
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

function post(op::LazyTextbookDiscretePost,
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
    tube⋂inv = Vector{ReachSet{LazySet{N}, N}}()
    source_invariant = HS.modes[source_loc_id].X
    dirs = op.options[:overapproximation]
    for reach_set in reach_tube
        R⋂I = Intersection(reach_set.X, source_invariant)
        if op.options[:check_invariant_intersection] && isempty(R⋂I)
            break
        end
        if !op.options[:lazy_R⋂I]
            R⋂I = overapproximate(R⋂I, dirs)
        end
        push!(tube⋂inv, ReachSet{LazySet{N}, N}(R⋂I,
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
    source_invariant = HS.modes[source_loc_id].X
    inv_isa_Hrep, inv_isa_H_polytope = get_Hrep_info(source_invariant)
    for trans in out_transitions(HS, source_loc_id)
        info("Considering transition: $trans")
        target_loc_id = target(HS, trans)
        target_loc = HS.modes[target(HS, trans)]
        target_invariant = target_loc.X
        trans_annot = HS.resetmaps[symbol(HS, trans)]
        guard = trans_annot.X
        assignment = trans_annot.A

        if inv_isa_Hrep
            guard_isa_Hrep, guard_isa_H_polytope = get_Hrep_info(guard)
        end

        # perform jumps
        post_jump = Vector{ReachSet{LazySet{N}, N}}()
        sizehint!(post_jump, length(tube⋂inv))
        for reach_set in tube⋂inv
            # check intersection with guard
            taken_intersection = false
            if inv_isa_Hrep && guard_isa_Hrep && op.options[:lazy_R⋂I]
                # combine the constraints of invariant and guard
                T = inv_isa_H_polytope || guard_isa_H_polytope ?
                    HPolytope :
                    HPolyhedron
                invariant_guard = T([constraints_list(source_invariant);
                    constraints_list(guard)])
                R⋂G = Intersection(reach_set.X.X, invariant_guard)
                taken_intersection = true
            end
            if !taken_intersection
                R⋂G = Intersection(reach_set.X, guard)
            end
            if isempty(R⋂G)
                continue
            end

            # apply assignment
            A⌜R⋂G⌟ = LinearMap(assignment, R⋂G)
            if !op.options[:lazy_R⋂G]
               A⌜R⋂G⌟ = overapproximate(A⌜R⋂G⌟, dirs)
            end

            # intersect with target invariant
            A⌜R⋂G⌟⋂I = Intersection(target_invariant, A⌜R⋂G⌟)

            # check if the final set is empty
            if isempty(A⌜R⋂G⌟⋂I)
                continue
            end

            # overapproximate final set once more
            if !op.options[:lazy_A⌜R⋂G⌟⋂I]
                res = overapproximate(A⌜R⋂G⌟⋂I, dirs)
            else
                res = A⌜R⋂G⌟⋂I
            end



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
