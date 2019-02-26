export ConcreteBlocksDiscretePost

"""
    ConcreteBlocksDiscretePost <: DiscretePost

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
struct ConcreteBlocksDiscretePost <: DiscretePost
    options::Options

    function ConcreteBlocksDiscretePost(𝑂::Options)
        𝑂copy = copy(𝑂)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:check_invariant_intersection], false)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:overapproximation], Hyperrectangle)
        return new(𝑂copy)
    end
end

# convenience constructor from pairs of symbols
ConcreteBlocksDiscretePost(𝑂::Pair{Symbol,<:Any}...) = ConcreteBlocksDiscretePost(Options(Dict{Symbol,Any}(𝑂)))

# default options for the LazyDiscretePost discrete post operator
ConcreteBlocksDiscretePost() = ConcreteBlocksDiscretePost(Options())

init(𝒫::ConcreteBlocksDiscretePost, 𝒮::AbstractSystem, 𝑂::Options) = init!(𝒫, 𝒮, copy(𝑂))

function init!(𝒫::ConcreteBlocksDiscretePost, 𝒮::AbstractSystem, 𝑂::Options)
    @assert isdefined(Main, :Polyhedra) "this algorithm needs the package " *
            "'Polyhedra' to be loaded"

    𝑂.dict[:n] = statedim(𝒮, 1)

    # solver-specific options (adds default values for unspecified options)
    𝑂out = validate_solver_options_and_add_default_values!(𝑂)

    # Input -> Output variable mapping
    𝑂out.dict[:inout_map] = inout_map_reach(𝑂out[:partition], 𝑂out[:blocks], 𝑂out[:n])
    return 𝑂out
end

function tube⋂inv!(𝒫::ConcreteBlocksDiscretePost,
                   reach_tube::Vector{<:ReachSet{<:LazySet{N}}},
                   invariant,
                   Rsets,
                   start_interval
                  ) where {N}
    # take intersection with source invariant
    dirs = 𝒫.options[:overapproximation]
    # TODO First check for empty intersection, which can be more efficient.
    #      However, we need to make sure that the emptiness check does not just
    #      compute the concrete intersection; otherwise, we would do the work
    #      twice. This is currently the case for 'Polyhedra' polytopes.

    # counts the number of sets R⋂I added to Rsets
    count = 0
    for reach_set in reach_tube

        R⋂I = intersection(invariant, reach_set.X)
        if isempty(R⋂I)
            break
        end
        #R⋂I = overapproximate(R⋂I, dirs)
        push!(Rsets, ReachSet{LazySet{N}, N}(R⋂I,
            reach_set.t_start + start_interval[1],
            reach_set.t_end + start_interval[2]))
        count = count + 1
    end

    return count
end

function post(𝒫::ConcreteBlocksDiscretePost,
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
        sizehint!(post_jump, count_Rsets)
        println(length(tube⋂inv))
        for reach_set in tube⋂inv[length(tube⋂inv) - count_Rsets + 1 : end]
            # check intersection with guard
            R⋂G = intersection(guard, reach_set.X)
            if isempty(R⋂G)
                continue
            end

            # apply assignment
            # TODO converting to HPolytope ?? handle automatically ??

            #A⌜R⋂G⌟ = LinearMap(assignment, oR)

            A⌜R⋂G⌟ = linear_map(assignment, R⋂G)

            # intersect with target invariant
            A⌜R⋂G⌟⋂I = intersection(A⌜R⋂G⌟, target_invariant)

            if isempty(A⌜R⋂G⌟⋂I)
                continue
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
