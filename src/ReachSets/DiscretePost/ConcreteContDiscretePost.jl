export ConcreteContDiscretePost
import Reachability.solve!
"""
    ConcreteContDiscretePost <: DiscretePost

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
struct ConcreteContDiscretePost <: DiscretePost
    options::Options

    function ConcreteContDiscretePost(𝑂::Options)
        𝑂copy = copy(𝑂)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:check_invariant_intersection], false)
        check_aliases_and_add_default_value!(𝑂.dict, 𝑂copy.dict, [:overapproximation], Hyperrectangle)
        return new(𝑂copy)
    end
end

# convenience constructor from pairs of symbols
ConcreteContDiscretePost(𝑂::Pair{Symbol,<:Any}...) = ConcreteContDiscretePost(Options(Dict{Symbol,Any}(𝑂)))

# default options for the LazyDiscretePost discrete post operator
ConcreteContDiscretePost() = ConcreteContDiscretePost(Options())

init(𝒫::ConcreteContDiscretePost, 𝒮::AbstractSystem, 𝑂::Options) = init!(𝒫, 𝒮, copy(𝑂))

function init!(𝒫::ConcreteContDiscretePost, 𝒮::AbstractSystem, 𝑂::Options)
    @assert isdefined(Main, :Polyhedra) "this algorithm needs the package " *
            "'Polyhedra' to be loaded"

    𝑂.dict[:n] = statedim(𝒮, 1)

    # solver-specific options (adds default values for unspecified options)
    𝑂out = validate_solver_options_and_add_default_values!(𝑂)

    return 𝑂out
end

function tube⋂inv!(𝒫::ConcreteContDiscretePost,
                   reach_tube::Vector{<:ReachSet{<:LazySet{N}}},
                   invariant,
                   Rsets,
                   low_dim_sets,
                   start_interval,
                   global_vars,
                   computed_vars) where {N}

    dirs = 𝒫.options[:overapproximation]

    # counts the number of sets R⋂I added to Rsets
    count = 0
    invariant_proj = HPolytope(constraints_list(LazySets.Approximations.project(invariant, computed_vars, LinearMap)))
    for reach_set in reach_tube
        R⋂I = intersection(reach_set.X, invariant_proj, true)
        if isempty(R⋂I)
            break
        end
        push!(Rsets, ReachSet{LazySet{N}, N}(Approximations.project(R⋂I, global_vars, dirs),
            reach_set.t_start + start_interval[1],
            reach_set.t_end + start_interval[2],
            reach_set.k))
        push!(low_dim_sets, ReachSet{LazySet{N}, N}(R⋂I,
            reach_set.t_start + start_interval[1],
            reach_set.t_end + start_interval[2],
            reach_set.k))
        count = count + 1
    end

    return count
end

function post(𝒫::ConcreteContDiscretePost,
              opC::ContinuousPost,
              HS::HybridSystem,
              waiting_list::Vector{Tuple{Int, ReachSet{LazySet{N}, N}, Int}},
              passed_list,
              source_loc_id,
              tube⋂inv,
              low_dim_sets,
              count_Rsets,
              jumps,
              options,
              options_copy,
              X0,
              computed_vars,
              global_vars) where {N}
    jumps += 1
    dirs = get_overapproximation_option(𝒫, options[:n])

    post_fix = Dict{Int, Vector{ReachSet{LazySet{N}, N}}}()
    steps = Vector{Int}()
    sizehint!(steps, count_Rsets)
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
        guard_proj = HPolytope(constraints_list(LazySets.Approximations.project(guard, computed_vars, LinearMap)))
        target_invariant_proj = HPolytope(constraints_list(LazySets.Approximations.project(target_invariant, computed_vars, LinearMap)))
        for reach_set in low_dim_sets
            # check intersection with guard

            R⋂G = intersection(reach_set.X,guard_proj,true)
            if isempty(R⋂G)
                continue
            end
            A⌜R⋂G⌟ = linear_map(assignment, R⋂G, computed_vars)
            A⌜R⋂G⌟⋂I = intersection(A⌜R⋂G⌟, target_invariant_proj, true)
            if isempty(A⌜R⋂G⌟⋂I)
                continue
            end
            push!(post_jump, ReachSet{LazySet{N}, N}(Approximations.project(A⌜R⋂G⌟⋂I, global_vars, dirs),
                                                     reach_set.t_start,
                                                     reach_set.t_end, reach_set.k))
            push!(steps, reach_set.k)
        end

        post_fix[length(post_fix)+1] = post_jump
    end

    unique!(steps)

    println("Start_high")
    println(length(steps))
    loc = HS.modes[source_loc_id]
    opc_options_copy = copy(opC.options.specified)
    opc_options_copy.dict[:vars] = 1:statedim(loc)
    opc_options_copy.dict[:steps] = steps
    println("Start_high")
    opC = BFFPSV18(opc_options_copy)

    h_reach_tube = solve!(IVP(loc, X0.X),
                        options_copy,
                        op=opC)

    h_Rsets = Vector{ReachSet{LazySet{N}, N}}()
    count = 1

    @inbounds for i = 1:length(h_reach_tube.Xk)
       reach_set = h_reach_tube.Xk[i]
       R⋂I = intersection(reach_set.X, loc.X, true)
       if isempty(R⋂I)
           break
       end
       push!(h_Rsets, ReachSet{LazySet{N}, N}(R⋂I,
           reach_set.t_start + X0.t_start,
           reach_set.t_end + X0.t_end,
           reach_set.k))

   end

    sizehint!(h_Rsets, length(steps))
    println(length(steps))
    println(length(h_reach_tube.Xk))



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
            R⋂G = intersection(reach_set.X, guard, true)
            if isempty(R⋂G)
                continue
            end

            #A⌜R⋂G⌟ = LinearMap(assignment, oR)
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
        end

        println("postprocess")

        postprocess(𝒫, HS, post_jump, options, waiting_list, passed_list,
            target_loc_id, jumps, post_fix[count])
        println("2postprocess")

        count = count + 1
    end

end
