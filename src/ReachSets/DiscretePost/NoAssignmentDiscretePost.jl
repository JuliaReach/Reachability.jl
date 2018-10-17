# ===========================================================================
# A discrete post operator for models without assignments and with polyhedral
# constraints (invariants and guards) only.
# ===========================================================================

import LazySets.use_precise_ρ

struct NoAssignmentDiscretePost <: DiscretePost
    options::Options
end

# default options for the NoAssignmentDiscretePost discrete post operator
function NoAssignmentDiscretePost()
    defaults = Options()
    setindex!(defaults, Hyperrectangle, :overapproximation)
    setindex!(defaults, false, :check_invariant_intersection)
    return NoAssignmentDiscretePost(defaults)
end

function init(op::NoAssignmentDiscretePost, HS::HybridSystem, options_input)
    options_input.dict[:n] = statedim(HS, 1)

    # solver-specific options (adds default values for unspecified options)
    options = validate_solver_options_and_add_default_values!(options_input)

    # Input -> Output variable mapping
    options.dict[:inout_map] =
        inout_map_reach(options[:partition], options[:blocks], options[:n])

    # set up operator-specific options
    @assert haskey(op.options.dict, :overapproximation)

    # TODO insert data into map
    op.options.dict[:set_map] = preprocess_model(op, HS)

    return options
end

function post(op::NoAssignmentDiscretePost,
              HS::HybridSystem,
              waiting_list::Vector{Tuple{Int, ReachSet{LazySet{N}, N}, Int}},
              passed_list,
              reach_tube::Vector{<:ReachSet{<:LazySet{N}}},
              Rsets,
              start_interval,
              source_loc_id,
              jumps,
              options
             ) where {N}
    # terminate if no jump should be taken anymore
    if jumps == options[:max_jumps]
        # take intersection with source invariant
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

        return nothing
    end

    jumps += 1
    dirs = get_overapproximation_option(op, options[:n])
    set_map = op.options[:set_map]
    source_invariant = HS.modes[source_loc_id].X
    for trans in out_transitions(HS, source_loc_id)
        info("Considering transition: $trans")
        target_loc_id = target(HS, trans)
        trans_id = symbol(HS, trans)
        set = set_map[source_loc_id][trans_id]

        # perform jumps
        post_jump = Vector{ReachSet{LazySet{N}, N}}()
        sizehint!(post_jump, length(reach_tube))
        for reach_set in reach_tube
            R⋂S = Intersection(reach_set.X, set)

            if isempty(R⋂S)
                if op.options[:check_invariant_intersection] &&
                        isempty(Intersection(reach_set.X, source_invariant))
                    break
                end
                continue
            end

            # overapproximate final set
            res = overapproximate(R⋂S, dirs)

            # store result
            push!(post_jump, ReachSet{LazySet{N}, N}(res,
                                                     reach_set.t_start,
                                                     reach_set.t_end))
        end

        postprocess(op, HS, post_jump, options, waiting_list, passed_list,
            target_loc_id, jumps)
    end
end

function preprocess_model(op::NoAssignmentDiscretePost, HS::HybridSystem)
    set_map = Vector{Dict{Int, LazySet}}(length(HS.modes))
    op.options.dict[:set_map] = set_map
    for source_loc_id in 1:length(HS.modes)
        inner_map = Dict{Int, LazySet}()
        set_map[source_loc_id] = inner_map
        loc = HS.modes[source_loc_id]
        source_invariant = loc.X
        s_inv_isa_Hrep, s_inv_isa_H_polytope = get_Hrep_info(source_invariant)
        for trans in out_transitions(HS, source_loc_id)
            trans_id = symbol(HS, trans)
            trans_annot = HS.resetmaps[trans_id]
            guard = trans_annot.X
            assignment = trans_annot.A
            target_loc_id = target(HS, trans)
            target_loc = HS.modes[target(HS, trans)]
            target_invariant = target_loc.X
            guard_isa_Hrep, guard_isa_H_polytope = get_Hrep_info(guard)
            t_inv_isa_Hrep, t_inv_isa_H_polytope =
                get_Hrep_info(target_invariant)
            if s_inv_isa_Hrep && guard_isa_Hrep && t_inv_isa_Hrep &&
                    assignment == eye(assignment)
                # simple model, unify constraints
                T = (s_inv_isa_H_polytope && guard_isa_H_polytope &&
                        t_inv_isa_H_polytope) ? HPolytope : HPolyhedron
                set = T([constraints_list(source_invariant);
                         constraints_list(guard);
                         constraints_list(target_invariant)]);
            else
                error("the model is not simple enough for the operator " *
                      "$(typeof(op))")
            end

            inner_map[trans_id] = set
        end
    end
    return set_map
end

# --- line search policies ---

# usually do not use line search
function use_precise_ρ(op::NoAssignmentDiscretePost,
                             cap::Intersection{N})::Bool where N<:Real
    return false
end

# use line search for the outermost level, which is a LinearMap
function use_precise_ρ(op::NoAssignmentDiscretePost,
                             cap::Intersection{N, <:LinearMap{N}}
                            )::Bool where N<:Real
    return true
end
