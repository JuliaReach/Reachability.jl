export ConcreteMixedSetDiscretePost

"""
    ConcreteMixedSetDiscretePost <: DiscretePost

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
struct ConcreteMixedSetDiscretePost <: DiscretePost
    options::Options

    function ConcreteMixedSetDiscretePost(𝑂::Options)
        𝑂copy = copy(𝑂)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:check_invariant_intersection], false)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:overapproximation], Hyperrectangle)
        return new(𝑂copy)
    end
end

# convenience constructor from pairs of symbols
ConcreteMixedSetDiscretePost(𝑂::Pair{Symbol,<:Any}...) = ConcreteMixedSetDiscretePost(Options(Dict{Symbol,Any}(𝑂)))

# default options for the LazyDiscretePost discrete post operator
ConcreteMixedSetDiscretePost() = ConcreteMixedSetDiscretePost(Options())

init(𝒫::ConcreteMixedSetDiscretePost, 𝒮::AbstractSystem, 𝑂::Options) = init!(𝒫, 𝒮, copy(𝑂))

function init!(𝒫::ConcreteMixedSetDiscretePost, 𝒮::AbstractSystem, 𝑂::Options)
    @assert isdefined(Main, :Polyhedra) "this algorithm needs the package " *
            "'Polyhedra' to be loaded"

    𝑂.dict[:n] = statedim(𝒮, 1)

    # solver-specific options (adds default values for unspecified options)
    𝑂out = validate_solver_options_and_add_default_values!(𝑂)

    return 𝑂out
end

function tube⋂inv!(𝒫::ConcreteMixedSetDiscretePost,
                   reach_tube::Vector{<:ReachSet{<:LazySet{N}}},
                   invariant,
                   Rsets,
                   low_dim_sets,
                   start_interval,
                   global_vars,
                   computed_vars
                  ) where {N}
    # take intersection with source invariant
    dirs = 𝒫.options[:overapproximation]

    # counts the number of sets R⋂I added to Rsets
    count = 0
    for reach_set in reach_tube
        # R⋂I = intersection(reach_set.X, invariant, true)
        # if isempty(R⋂I)
        #     break
        # end
        #R⋂I = overapproximate(R⋂I, dirs)
        push!(Rsets, ReachSet{LazySet{N}, N}(LazySets.Approximations.project(reach_set.X,global_vars,LinearMap),
            reach_set.t_start + start_interval[1],
            reach_set.t_end + start_interval[2],
            reach_set.k))
        if (dim(reach_set.X) != length(computed_vars))
            push!(low_dim_sets, ReachSet{LazySet{N}, N}(reach_set.X,
                reach_set.t_start + start_interval[1],
                reach_set.t_end + start_interval[2],
                reach_set.k))
        end
        count = count + 1
    end

    return count
end

function post(𝒫::ConcreteMixedSetDiscretePost,
                  HS::HybridSystem,
                  waiting_list::Vector{Tuple{Int, ReachSet{LazySet{N}, N}, Int}},
                  passed_list,
                  source_loc_id,
                  tube⋂inv,
                  low_dim_sets,
                  count_Rsets,
                  jumps,
                  options,
                  global_vars) where {N}
    jumps += 1
    dirs = get_overapproximation_option(𝒫, options[:n])
    for trans in out_transitions(HS, source_loc_id)
        println("trans")
        info("Considering transition: $trans")
        target_loc_id = target(HS, trans)
        target_loc = HS.modes[target(HS, trans)]
        target_invariant = target_loc.X
        trans_annot = HS.resetmaps[symbol(HS, trans)]
        guard = trans_annot.X
        assignment = trans_annot.A

        # perform jumps
        post_jump = Vector{ReachSet{LazySet{N}, N}}()
        post_jump_low = Vector{ReachSet{LazySet{N}, N}}()
        sizehint!(post_jump, count_Rsets)
        sizehint!(post_jump_low, count_Rsets)
        for reach_set in low_dim_sets
            # check intersection with guard
            R⋂G = intersection(reach_set.X, guard, true)
            if isempty(R⋂G)
                continue
            end

            # apply assignment

            A⌜R⋂G⌟ = linear_map(assignment, R⋂G)

            # intersect with target invariant
            A⌜R⋂G⌟⋂I = intersection(A⌜R⋂G⌟, target_invariant, true)

            if isempty(A⌜R⋂G⌟⋂I)
                continue
            end

            # store result
            push!(post_jump, ReachSet{LazySet{N}, N}(A⌜R⋂G⌟⋂I,
                                                     reach_set.t_start,
                                                     reach_set.t_end, reach_set.k))
            push!(post_jump_low, ReachSet{LazySet{N}, N}(LazySets.Approximations.project(A⌜R⋂G⌟⋂I, global_vars, LinearMap),
                                                  reach_set.t_start,
                                                  reach_set.t_end, reach_set.k))

        end

        println(length(post_jump))

        postprocess(𝒫, HS, post_jump, options, waiting_list, passed_list,
            target_loc_id, jumps, post_jump_low)
    end
end
