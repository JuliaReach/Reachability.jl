export LazyContDiscretePost

import LazySets.use_precise_ρ
import Reachability.solve!

"""
    LazyContDiscretePost <: DiscretePost

Textbook implementation of a discrete post operator, but with lazy intersections in low-dimensional systems first.

### Fields

- `options` -- an `Options` structure that holds the algorithm-specific options

### Algorithm

The algorithm is based on [Flowpipe-Guard Intersection for Reachability
Computations with Support Functions](http://spaceex.imag.fr/sites/default/files/frehser_adhs2012.pdf).
"""
struct LazyContDiscretePost <: DiscretePost
    options::Options

    function LazyContDiscretePost(𝑂::Options)
        𝑂copy = copy(𝑂)
        # TODO: pass 𝑂 directly?
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:check_invariant_intersection], true)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:overapproximation], Hyperrectangle)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:lazy_R⋂I], true)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:lazy_R⋂G], false)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:lazy_A⌜R⋂G⌟⋂I], false)
        return new(𝑂copy)
    end
end

# convenience constructor from pairs of symbols
LazyContDiscretePost(𝑂::Pair{Symbol,<:Any}...) = LazyContDiscretePost(Options(Dict{Symbol,Any}(𝑂)))

# default options for the LazyContDiscretePost discrete post operator
LazyContDiscretePost() = LazyContDiscretePost(Options())

init(𝒫::LazyContDiscretePost, 𝒮::AbstractSystem, 𝑂::Options) = init!(𝒫, 𝒮, copy(𝑂))

# TODO: use 𝑂 only?
function init!(𝒫::LazyContDiscretePost, 𝒮::AbstractSystem, 𝑂::Options)
    𝑂[:n] = statedim(𝒮, 1)

    # solver-specific options (adds default values for unspecified options)
    𝑂out = validate_solver_options_and_add_default_values!(𝑂)

    # Input -> Output variable mapping
    𝑂out[:inout_map] = inout_map_reach(𝑂out[:partition], 𝑂out[:blocks], 𝑂out[:n])

    return 𝑂out
end

function tube⋂inv!(𝒫::LazyContDiscretePost,
                   reach_tube::Vector{<:ReachSet{<:LazySet{N}}},
                   invariant,
                   Rsets,
                   start_interval,
                   options,
                   nonzero_vars
                  ) where {N}

    # TODO dirs = get_overapproximation_option(op, dim(invariant)) ?
    dirs = 𝒫.options[:overapproximation]

    # counts the number of sets R⋂I added to Rsets
    count = 0
    @inbounds for i = 1:length(reach_tube)
        invariant_proj = LazySets.Approximations.project(invariant, nonzero_vars, LinearMap)
        #Rimage_proj = LazySets.Approximations.overapproximate(reach_tube[i].X)

        R⋂I = Intersection(reach_tube[i].X, invariant_proj)

        if 𝒫.options[:check_invariant_intersection] && isempty(R⋂I)
            break
        end
        if !𝒫.options[:lazy_R⋂I]
            R⋂I = overapproximate(R⋂I, dirs)
        end

        push!(Rsets, ReachSet{LazySet{N}, N}(R⋂I,
            reach_tube[i].t_start + start_interval[1],
            reach_tube[i].t_end + start_interval[2],
            reach_tube[i].k))
        count = count + 1
    end

    return count
end

function post(𝒫::LazyContDiscretePost,
              HS::HybridSystem,
              waiting_list::Vector{Tuple{Int, ReachSet{LazySet{N}, N}, Int}},
              passed_list,
              source_loc_id,
              tube⋂inv,
              count_Rsets,
              jumps,
              options,
              nonzero_vars,
              opC,
              X0,
              time_horizon) where {N}
    # TODO? dirs = 𝒫.options[:overapproximation]
    dirs = get_overapproximation_option(𝒫, options[:n])
    source_invariant = HS.modes[source_loc_id].X
    inv_isa_Hrep, inv_isa_H_polytope = get_Hrep_info(source_invariant)
    steps = Vector{Int}()

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
        combine_constraints = inv_isa_Hrep && guard_isa_Hrep && 𝒫.options[:lazy_R⋂I]
        if combine_constraints # combine the constraints of invariant and guard
            T = inv_isa_H_polytope || guard_isa_H_polytope ? HPolytope : HPolyhedron
            # TODO: remove redundant constraints => use intersection(..)
            invariant_guard = T([constraints_list(source_invariant);
                                 constraints_list(guard)])

        end

        guard_proj = combine_constraints ? LazySets.Approximations.project(invariant_guard, nonzero_vars, LinearMap) : LazySets.Approximations.project(guard, nonzero_vars, LinearMap)

        # perform jumps
        post_jump = Vector{ReachSet{LazySet{N}, N}}()
        to_high_res = Vector{ReachSet{LazySet{N}, N}}()
        sizehint!(post_jump, count_Rsets)
        sizehint!(to_high_res, count_Rsets)
        #tube⋂inv[length(tube⋂inv) - count_Rsets + 1 : end]
        for reach_set in tube⋂inv[length(tube⋂inv) - count_Rsets + 1 : end]
            taken_intersection = false
            if combine_constraints
                    R⋂G = Intersection(reach_set.X.X, guard_proj)
                    taken_intersection = true
            end
            if !taken_intersection
                    R⋂G = Intersection(reach_set.X, guard_proj)
            end
            if isempty(R⋂G)
                continue
            end

            # apply assignment
            #R⋂G = overapproximate(R⋂G, dirs)
            push!(steps, reach_set.k)
        end
    end

    loc = HS.modes[source_loc_id]
    options_copy = copy(options)
    options_copy_low = copy(options)
    options_copy.dict[:T] = time_horizon - X0.t_start
    options_copy.dict[:project_reachset] = false
    delete!(options_copy.dict, :inout_map)
    #delete_N = !haskey(options_input, :N)
    if true # TODO add more conditions or fix option clashes in general
        delete!(options_copy.dict, :N)
    end
    if haskey(options_copy, :block_types) &&
            options_copy.dict[:block_types] == nothing
        delete!(options_copy.dict, :block_types)
    end
    if haskey(options_copy, :blocks)
        delete!(options_copy.dict, :blocks)
    end
    options_copy.dict[:vars] = 1:statedim(loc)
    options_copy.dict[:steps] = steps
    h_reach_tube = solve!(ContinuousSystem(loc.A, X0.X, loc.U),
                        options_copy,
                        op=opC,
                        invariant=source_invariant)
    h_Rsets = Vector{ReachSet{LazySet{N}, N}}()
    count = 0

    @inbounds for i = 1:length(h_reach_tube.Xk)
        reach_set = h_reach_tube.Xk[i]
        R⋂I = intersection(loc.X, reach_set.X)
        if isempty(R⋂I)
            break
        end
        #R⋂I = overapproximate(R⋂I, dirs)
        push!(h_Rsets, ReachSet{LazySet{N}, N}(R⋂I,
            reach_set.t_start + X0.t_start,
            reach_set.t_end + X0.t_end,
            reach_set.k))
        count = count + 1
    end
    jumps += 1
    dirs = get_overapproximation_option(𝒫, options[:n])
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
        sizehint!(post_jump, count)


        for r_i in 1:length(h_Rsets)
            reach_set = h_Rsets[r_i]
            #println(reach_set)
            # check intersection with guard
            R⋂G = intersection(guard, reach_set.X)
            if isempty(R⋂G)
                continue
            end

            #A⌜R⋂G⌟ = LinearMap(assignment, oR)
            A⌜R⋂G⌟ = linear_map(assignment, R⋂G)

            # intersect with target invariant
            A⌜R⋂G⌟⋂I = intersection(A⌜R⋂G⌟, target_invariant)

            if isempty(A⌜R⋂G⌟⋂I)
                continue
            end

            # store result
            push!(post_jump, ReachSet{LazySet{N}, N}(overapproximate(A⌜R⋂G⌟⋂I),
                                                     reach_set.t_start,
                                                     reach_set.t_end, reach_set.k))
        end
        println(length(post_jump))
        postprocess(𝒫, HS, post_jump, options, waiting_list, passed_list,
            target_loc_id, jumps)
    end
end

function get_Hrep_info(set::LazySet)
    return (false, false)
end

function get_Hrep_info(set::HPolytope)
    return (true, true)
end

function get_Hrep_info(set::HPolyhedron)
    return (true, false)
end

# --- line search policies ---

# usually do not use line search
function use_precise_ρ(𝒫::LazyContDiscretePost,
                       cap::Intersection{N})::Bool where N<:Real
    return false
end

# use line search for the outermost level, which is a LinearMap
function use_precise_ρ(𝒫::LazyContDiscretePost,
                       cap::Intersection{N, <:LinearMap{N}}
                       )::Bool where N<:Real
    return true
end
