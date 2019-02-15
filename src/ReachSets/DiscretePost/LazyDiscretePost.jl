export LazyDiscretePost,
       ApproximatingDiscretePost

import LazySets.use_precise_ρ

"""
    LazyDiscretePost <: DiscretePost

Textbook implementation of a discrete post operator, but with lazy intersections.

### Fields

- `options` -- an `Options` structure that holds the algorithm-specific options

### Algorithm

The algorithm is based on [Flowpipe-Guard Intersection for Reachability
Computations with Support Functions](http://spaceex.imag.fr/sites/default/files/frehser_adhs2012.pdf).
"""
struct LazyDiscretePost <: DiscretePost
    options::Options

    function LazyDiscretePost(𝑂::Options)
        𝑂copy = copy(𝑂)
        # TODO: pass 𝑂 directly?
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:check_invariant_intersection], false)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:overapproximation], Hyperrectangle)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:lazy_R⋂I], false)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:lazy_R⋂G], true)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:lazy_A⌜R⋂G⌟], true)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:lazy_A⌜R⋂G⌟⋂I], true)
        return new(𝑂copy)
    end
end

# convenience constructor from pairs of symbols
LazyDiscretePost(𝑂::Pair{Symbol,<:Any}...) = LazyDiscretePost(Options(Dict{Symbol,Any}(𝑂)))

# default options for the LazyDiscretePost discrete post operator
LazyDiscretePost() = LazyDiscretePost(Options())

"""
    ApproximatingDiscretePost()

Textbook implementation of a discrete post operator, but with lazy intersections
followed by an overapproximation. This is a particular case of the
`LazyDiscretePost`.
"""
function ApproximatingDiscretePost()
    return LazyDiscretePost(:check_invariant_intersection=>false,
                            :overapproximation=>Hyperrectangle,
                            :lazy_R⋂I=>false,
                            :lazy_R⋂G=>false,
                            :lazy_A⌜R⋂G⌟⋂I=>false)
end

function ApproximatingDiscretePost(𝑂::Options)
    𝑂_default = Options(:lazy_R⋂I=>false,
                        :lazy_R⋂G=>false,
                        :lazy_A⌜R⋂G⌟⋂I=>false)
    merge!(𝑂_default, 𝑂)
    LazyDiscretePost(𝑂_default)
end

init(𝒫::LazyDiscretePost, 𝒮::AbstractSystem, 𝑂::Options) = init!(𝒫, 𝒮, copy(𝑂))

# TODO: use 𝑂 only?
function init!(𝒫::LazyDiscretePost, 𝒮::AbstractSystem, 𝑂::Options)
    𝑂[:n] = statedim(𝒮, 1)

    # solver-specific options (adds default values for unspecified options)
    𝑂out = validate_solver_options_and_add_default_values!(𝑂)

    # Input -> Output variable mapping
    𝑂out[:inout_map] = inout_map_reach(𝑂out[:partition], 𝑂out[:blocks], 𝑂out[:n])

    return 𝑂out
end

function tube⋂inv!(𝒫::LazyDiscretePost,
                   reach_tube::Vector{<:ReachSet{<:LazySet{N}}},
                   invariant,
                   Rsets,
                   start_interval
                  ) where {N}

    # TODO dirs = get_overapproximation_option(op, dim(invariant)) ?
    dirs = 𝒫.options[:overapproximation]

    # counts the number of sets R⋂I added to Rsets
    count = 0
    @inbounds for reach_set in reach_tube
        R⋂I = Intersection(reach_set.X, invariant)
        if 𝒫.options[:check_invariant_intersection] && isempty(R⋂I)
            break
        end
        if !𝒫.options[:lazy_R⋂I]
            R⋂I = overapproximate(R⋂I, dirs)
        end
        push!(Rsets, ReachSet{LazySet{N}, N}(R⋂I,
            reach_set.t_start + start_interval[1],
            reach_set.t_end + start_interval[2]))
        count = count + 1
    end

    return count
end

function post(𝒫::LazyDiscretePost,
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
    # TODO? oa = 𝒫.options[:overapproximation]
    oa = get_overapproximation_option(𝒫, options[:n])
    source_invariant = HS.modes[source_loc_id].X
    inv_isa_Hrep, inv_isa_H_polytope = get_Hrep_info(source_invariant)

    for trans in out_transitions(HS, source_loc_id)
        info("Considering transition: $trans")
        target_loc_id = target(HS, trans)
        target_loc = HS.modes[target(HS, trans)]
        target_invariant = target_loc.X
        constrained_map = resetmap(HS, trans)
        guard = stateset(constrained_map)

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

        # perform jumps
        post_jump = Vector{ReachSet{LazySet{N}, N}}()
        sizehint!(post_jump, count_Rsets)
        for reach_set in tube⋂inv[length(tube⋂inv) - count_Rsets + 1 : end]
            # check intersection with guard
            if combine_constraints
                R⋂G = Intersection(reach_set.X.X, invariant_guard)
            else
                R⋂G = Intersection(reach_set.X, guard)
            end
            if isempty(R⋂G)
                continue
            end
            if !𝒫.options[:lazy_R⋂G]
                R⋂G = overapproximate(R⋂G, oa)
            end

            # apply assignment
            A⌜R⋂G⌟ = apply_assignment(𝒫, constrained_map, R⋂G)
            if !𝒫.options[:lazy_A⌜R⋂G⌟]
                A⌜R⋂G⌟ = overapproximate(A⌜R⋂G⌟, oa)
            end

            # intersect with target invariant
            A⌜R⋂G⌟⋂I = Intersection(target_invariant, A⌜R⋂G⌟)
            if isempty(A⌜R⋂G⌟⋂I)
                continue
            end
            if !𝒫.options[:lazy_A⌜R⋂G⌟⋂I]
                A⌜R⋂G⌟⋂I = overapproximate(A⌜R⋂G⌟⋂I, oa)
            end

            # store result
            push!(post_jump, ReachSet{LazySet{N}, N}(A⌜R⋂G⌟⋂I,
                                                     reach_set.t_start,
                                                     reach_set.t_end))
        end

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

# --- line-search policies ---

# usually do not use line search
function use_precise_ρ(𝒫::LazyDiscretePost,
                       cap::Intersection{N})::Bool where N<:Real
    return false
end

# use line search for the outermost level, which is a LinearMap
function use_precise_ρ(𝒫::LazyDiscretePost,
                       cap::Intersection{N, <:LinearMap{N}}
                       )::Bool where N<:Real
    return true
end
