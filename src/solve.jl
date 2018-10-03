using LazySets.Approximations.box_approximation
export AbstractSolution,
       ReachSolution,
       CheckSolution,
       solve,
       project

"""
    AbstractSolution

Abstract type representing the solution of a rechability problem.
"""
abstract type AbstractSolution end

"""
    ReachSolution{S<:LazySet} <: AbstractSolution

Type that wraps a the solution of a reachability problem as a sequence of
lazy sets, and a dictionary of options.

### Fields

- `Xk`       -- the list of reachable states
- `options`  -- the dictionary of options

### Notes

If the solution has been projected in 2D, the sequence `Xk` is an array
of polygons in constraint representation. In high-dimensions this is a sequence
of cartesian product arrays of low-dimensional sets.
"""
struct ReachSolution{S<:LazySet} <: AbstractSolution
  Xk::Vector{S}
  options::Options
end
# constructor with no options
ReachSolution(Xk::Vector{S}) where {S<:LazySet} =
    ReachSolution{S}(Xk, Options())

"""
    CheckSolution

Type that wraps a the solution of a property checking problem, which is just the
answer if the property is satisfied.

### Fields

- `satisfied` -- is the property satisfied?
- `violation` -- step at which the property is violated (-1 otherwise)
- `options`   -- the dictionary of options
"""
struct CheckSolution <: AbstractSolution
  satisfied::Bool
  violation::Int
  options::Options
end
# constructor with no options
CheckSolution(satisfied::Bool, violation::Int) =
    CheckSolution(satisfied, violation, Options())

function default_algorithm(system::InitialValueProblem)
    algorithm = ""
    s = system.s
    if s isa LinearContinuousSystem ||
       s isa LinearControlContinuousSystem ||
       s isa ConstrainedLinearContinuousSystem ||
       s isa ConstrainedLinearControlContinuousSystem ||
       s isa LinearDiscreteSystem ||
       s isa LinearControlDiscreteSystem ||
       s isa ConstrainedLinearDiscreteSystem ||
       s isa ConstrainedLinearControlDiscreteSystem
       
        algorithm = "BFFPSV18"
    else
        error("no default reachability algorithm available for system of " *
              "type $(typeof(system))")
    end
    return algorithm
end

"""
    solve(system, options)  or  solve(system, :key1 => val1, [...], keyK => valK)

Solves a reachability problem s.t. the given options.
If some options are not defined, we may fall back to default values.

### Input

- `system`    -- a (discrete or continuoues) system specification
- `options`   -- algorithm options for solving the problem
- `algorithm` -- (optional, default: dispatched on the system's type) the
                 reachability algorithm for the computation

### Output

A solution object whose content depends on the input options.

### Notes

To see all available input options, see
`keys(Reachability.available_keywords.dict)`.
"""
function solve(system::InitialValueProblem,
               options::Options;
               algorithm::String=default_algorithm(system))
    solve!(system, Options(copy(options.dict)), algorithm=algorithm)
end

solve(system::AbstractSystem, options::Pair{Symbol,<:Any}...) =
    solve(system, Options(Dict{Symbol,Any}(options)))

function solve!(system::InitialValueProblem{<:SYS},
                options::Options;
                algorithm::String=default_algorithm(system)
               )::AbstractSolution where {SYS<:Union{AbstractContinuousSystem,
                                                     AbstractDiscreteSystem}}
    if algorithm == "BFFPSV18"
        options = init_BFFPSV18!(system, options)

        # coordinate transformation
        options[:transformation_matrix] = nothing
        if options[:coordinate_transformation] != ""
            info("Transformation...")
            tic()
            (system, transformation_matrix) =
                transform(system, options[:coordinate_transformation])
            tocc()
            options[:transformation_matrix] = transformation_matrix
        end

        # convert matrix
        system = matrix_conversion(system, options)

        if options[:mode] == "reach"
            info("Reachable States Computation...")
            tic()
            Rsets = reach(system, options)
            info("- Total")
            tocc()

            # Projection
            if options[:project_reachset] || options[:projection_matrix] != nothing
                info("Projection...")
                tic()
                RsetsProj = project(Rsets, options)
                tocc()
                return ReachSolution(RsetsProj, options)
            end

            return ReachSolution(Rsets, options)

        elseif options[:mode] == "check"

            # Input -> Output variable mapping in property
            options.dict[:property] = inout_map_property(options[:property],
                options[:partition], options[:blocks], options[:n])

            # =================
            # Property checking
            # =================
            info("Property Checking...")
            tic()
            answer = check_property(system, options)
            info("- Total")
            tocc()

            if answer == 0
                info("The property is satisfied!")
                return CheckSolution(true, -1, options)
            else
                info("The property may be violated at index $answer," *
                    " (time point $(answer * options[:δ]))!")
                return CheckSolution(false, answer, options)
            end
        else
            error("unsupported mode $(options[:mode])")
        end # mode
    else
        error("unsupported algorithm $algorithm")
    end # algorithm
end

"""
    project(Rsets, options)

Projects a sequence of sets according to the settings defined in the options.

### Input

- `Rsets`   -- solution of a reachability problem
- `options` -- options structure

### Notes

A projection matrix can be given in the options structure, or passed as a
dictionary entry.
"""
function project(Rsets::Vector{<:LazySet}, options::Options)
    plot_vars = copy(options[:plot_vars])
    for i in 1:length(plot_vars)
        if plot_vars[i] != 0
            plot_vars[i] = options[:inout_map][plot_vars[i]]
        end
    end
    reduced_n = sum(x -> x != 0, options[:inout_map])
    output_function = !options[:project_reachset]
    RsetsProj = project_reach(plot_vars,
                              reduced_n,
                              options[:δ],
                              Rsets,
                              options[:algorithm],
                              ε=options[:ε_proj],
                              set_type=options[:set_type_proj],
                              transformation_matrix=options[:transformation_matrix],
                              projection_matrix=options[:projection_matrix],
                              output_function=output_function
                             )
end

project(reach_sol::AbstractSolution) = project(reach_sol.Xk, reach_sol.options)

project(Rsets::Vector{<:LazySet}, options::Pair{Symbol,<:Any}...) =
    project(Rsets, Options(Dict{Symbol,Any}(options)))

"""
    solve(HS::HybridSystem,
          X0::LazySet,
          options::Options)::AbstractSolution

Interface to reachability algorithms for a hybrid system PWA dynamics.

### Input

- `HS`      -- hybrid system
- `X0`      -- initial set
- `options` -- options for solving the problem

### Notes

The current implementation requires that you have loaded the `Polyhedra` library,
because some concrete operations between polytopes are used.

Currently, the following simplifying assumptions are made:

- the starting state (for which `X0` is given) correspond to the first location
  in the automaton
- source invariants, target invariants and guards are polytopes in constraint
  representation

### Algorithm

The algorithm is based on [Flowpipe-Guard Intersection for Reachability Computations
with Support Functions](http://spaceex.imag.fr/sites/default/files/frehser_adhs2012.pdf).
"""
function solve(HS::HybridSystem,
               X0::LazySet,
               options::Options)::AbstractSolution
    @assert isdefined(Main, :Polyhedra) "this algorithm needs the package " *
            "'Polyhedra' to be loaded"

    # set options
    options.dict[:project_reachset] = false

    max_jumps = options[:max_jumps]

    #TODO get initial location. For now we assume that it is the first location
    # waiting_list entries:
    # - (discrete) location
    # - (set of) continuous states
    # - number of previous jumps
    waiting_list::Vector{Tuple{Int, LazySet, Int}} = [(1, X0, 0)]
    rset = []
    while (!isempty(waiting_list))
        cur_loc_id, X0, jumps = pop!(waiting_list)
        cur_loc = HS.modes[cur_loc_id]
        S = ContinuousSystem(cur_loc.A, X0, cur_loc.U)

        # compute reach tubes
        Rsets = solve(S, options)

        cur_invariant = cur_loc.X

        # TODO temporary conversion to HPolytope
        @assert cur_invariant isa HalfSpace
        cur_invariant = HPolytope([cur_invariant])

        # take intersection with source invariant
        intersectedRset =
            intersect_reach_tubes_invariant(Rsets.Xk, cur_invariant)
        push!(rset, intersectedRset)

        if jumps == max_jumps
            continue
        end
        for trans in out_transitions(HS, cur_loc_id)
            info("Considering transition: $trans")
            target_loc_id = target(HS, trans)
            target_loc = HS.modes[target(HS, trans)]
            target_invariant = target_loc.X
            guard = HS.resetmaps[target_loc_id].X
            assignment = HS.resetmaps[target_loc_id].A

            # TODO temporary conversion to HPolytope
            @assert target_invariant isa HalfSpace
            target_invariant = HPolytope([target_invariant])
            @assert guard isa HPolytope

            # check intersection with guard
            rsetIntersMinus = [intersection(guard, VPolytope(vertices_list(hi))) for hi in intersectedRset]
            filter!(!isempty, rsetIntersMinus)
            if (isempty(rsetIntersMinus))
                continue
            end

            # apply assignment
            rsetIntersMinus = [linear_map(assignment, ri) for ri in rsetIntersMinus]

            # check intersection with target invariant
            info("Intersection with I\^+")
            rsetIntersPlus = [intersection(target_invariant, hi) for hi in rsetIntersMinus]
            filter!(!isempty, rsetIntersPlus)

            # push new sets after jump
            rsetHull = ConvexHullArray(rsetIntersPlus)
            for rh in rsetIntersPlus
                push!(waiting_list, (target(HS, trans), rh, jumps + 1))
            end

            # TODO check intersection with forbidden states
            # push!(waiting_list,(target(HS, trans), rsetHull))
        end
    end
    return ReachSolution(vcat(rset...), options)
end

function intersect_reach_tubes_invariant(reach_tubes, invariant)
    # TODO First check for empty intersection, which can be more efficient.
    #      However, we need to make sure that the emptiness check does not just
    #      compute the concrete intersection; otherwise, we would do the work
    #      twice. This is currently the case for 'Polyhedra' polytopes.
    # TODO pass numeric type correctly
    intersections = Vector{HPolytope{Float64}}()
    for rt in reach_tubes
        # TODO temporary workaround for 1D sets
        if dim(rt) == 1
            reach_tube = VPolytope(vertices_list(
                Approximations.overapproximate(rt, LazySets.Interval)))
        # TODO offer a lazy intersection here
        # TODO offer more options instead of taking the VPolytope intersection
        elseif rt isa CartesianProductArray
            reach_tube = VPolytope(vertices_list(rt))
        else
            error("unsupported set type for reach tube: $(typeof(rt))")
        end
        intersect = intersection(invariant, reach_tube)
        if isempty(intersect)
            break
        end
        push!(intersections, intersect)
    end
    return intersections
end

# ===========================================================================
# Bogomolov, Forets, Frehse, Podelski, Schilling, Viry 2018
# ===========================================================================
function init_BFFPSV18!(system, options_input)
    # state dimension for (purely continuous or purely discrete systems)
    options_input.dict[:n] = statedim(system)

    # solver-specific options (adds default values for unspecified options)
    options = validate_solver_options_and_add_default_values!(options_input)

    # Input -> Output variable mapping
    options.dict[:inout_map] = inout_map_reach(options[:partition], options[:blocks], options[:n])

    if options[:project_reachset]
        options[:output_function] = nothing
    else
        options[:output_function] = options[:projection_matrix]
    end

    return options
end
