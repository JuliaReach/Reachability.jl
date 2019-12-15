export ConcreteDiscretePost

"""
    ConcreteDiscretePost <: DiscretePost

Textbook implementation of a discrete post operator, using concrete polyhedra
intersections.

### Fields

- `options` -- an `Options` structure that holds the algorithm-specific options

### Notes

This operator requires that the `Polyhedra` library is loaded,
because some concrete operations between polytopes are used.

Currently, we assume that source invariants, target invariants and guards are
polytopes in constraint representation resp. half-spaces.

### Algorithm

The algorithm is based on [Flowpipe-Guard Intersection for Reachability
Computations with Support Functions](http://spaceex.imag.fr/sites/default/files/frehser_adhs2012.pdf).
"""
struct ConcreteDiscretePost <: DiscretePost
    options::Options

    function ConcreteDiscretePost(𝑂::Options)
        𝑂copy = copy(𝑂)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:check_invariant_intersection], false)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:overapproximation], Hyperrectangle)
        return new(𝑂copy)
    end
end

# convenience constructor from pairs of symbols
ConcreteDiscretePost(𝑂::Pair{Symbol,<:Any}...) = ConcreteDiscretePost(Options(Dict{Symbol,Any}(𝑂)))

# default options for the LazyDiscretePost discrete post operator
ConcreteDiscretePost() = ConcreteDiscretePost(Options())

init(𝒫::ConcreteDiscretePost, 𝒮::AbstractSystem, 𝑂::Options) = init!(𝒫, 𝒮, copy(𝑂))

function init!(𝒫::ConcreteDiscretePost, 𝒮::AbstractSystem, 𝑂::Options)
    @assert isdefined(Main, :Polyhedra) "this algorithm needs the package " *
            "'Polyhedra' to be loaded"

    𝑂.dict[:n] = statedim(𝒮, 1)

    # solver-specific options (adds default values for unspecified options)
    𝑂out = validate_solver_options_and_add_default_values!(𝑂)

    return 𝑂out
end

function tube⋂inv(𝒫::ConcreteDiscretePost,
                  reach_tube::Vector{<:AbstractReachSet{<:LazySet{N}}},
                  invariant,
                  start_interval
                 ) where {N}
    # take intersection with source invariant

    Rsets = Vector{AbstractReachSet{<:LazySet{N}}}()
    for reach_set in reach_tube
        rs = set(reach_set)
        if rs isa CartesianProductArray && length(array(rs)) == 1
            # TODO workaround for lazy X0
            rs = overapproximate(rs, BoxDirections)
        end
        if 𝒫.options[:check_invariant_intersection] && isdisjoint(rs, invariant)
            break
        end
        push!(Rsets,
              substitute(reach_set, set=intersection(rs, invariant),
                         time_start=time_start(reach_set) + start_interval[1],
                         time_end=time_end(reach_set) + start_interval[2]))
    end

    return Rsets
end

function post(𝒫::ConcreteDiscretePost,
              HS::HybridSystem,
              waiting_list::Vector{Tuple{Int, <:AbstractReachSet{LazySet{N}}, Int}},
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
        constrained_map = resetmap(HS, trans)
        guard = stateset(constrained_map)

        # perform jumps
        post_jump = Vector{ReachSet{LazySet{N}}}()
        sizehint!(post_jump, length(tube⋂inv))
        for reach_set in tube⋂inv
            # check intersection with guard
            R⋂G = intersection(guard, set(reach_set))
            if isempty(R⋂G)
                continue
            end

            # apply assignment
            A⌜R⋂G⌟ = apply_assignment(𝒫, constrained_map, R⋂G)

            # intersect with target invariant
            A⌜R⋂G⌟⋂I = intersection(target_invariant, A⌜R⋂G⌟)
            if isempty(A⌜R⋂G⌟⋂I)
                continue
            end

            # store result
            push!(post_jump, ReachSet{LazySet{N}}(A⌜R⋂G⌟⋂I,
                                                     time_start(reach_set),
                                                     time_end(reach_set)))
        end

        postprocess(𝒫, HS, post_jump, options, waiting_list, passed_list,
            target_loc_id, jumps)
    end
end

# --- handling assignments ---

function apply_assignment(𝒫::ConcreteDiscretePost,
                          constrained_map::ConstrainedLinearMap,
                          R⋂G::LazySet;
                          kwargs...)
    return linear_map(constrained_map.A, R⋂G)
end
