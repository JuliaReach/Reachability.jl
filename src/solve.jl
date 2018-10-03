export solve,
       project

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
            else
                RsetsProj = Rsets
            end

            return ReachSolution(RsetsProj, options)

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
function project(Rsets::Vector{<:ReachSet}, options::Options)
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

project(Rsets::Vector{<:ReachSet}, options::Pair{Symbol,<:Any}...) =
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

The current implementation requires that you have loaded the `Polyhedra`
library, because some concrete operations between polytopes are used.

Currently, the following simplifying assumptions are made:

- the starting state (for which `X0` is given) correspond to the first location
  in the automaton
- source invariants, target invariants and guards are polytopes in constraint
  representation

### Algorithm

The algorithm is based on [Flowpipe-Guard Intersection for Reachability
Computations with Support Functions](
http://spaceex.imag.fr/sites/default/files/frehser_adhs2012.pdf).
"""
function solve(HS::HybridSystem,
               X0::LazySet{N},
               options::Options)::AbstractSolution where {N}
    @assert isdefined(Main, :Polyhedra) "this algorithm needs the package " *
            "'Polyhedra' to be loaded"

    # set options
    options.dict[:project_reachset] = false

    max_jumps = options[:max_jumps]

    time_horizon = options[:T]

    #TODO get initial location. For now we assume that it is the first location
    # waiting_list entries:
    # - (discrete) location
    # - (set of) continuous-time reach sets
    # - number of previous jumps
    waiting_list::Vector{Tuple{Int, ReachSet{LazySet{N}, N}, Int}} =
        [(1, ReachSet{LazySet{N}, N}(X0, zero(N), zero(N)), 0)]
    rset = Vector{ReachSet{LazySet{N}, N}}()
    while (!isempty(waiting_list))
        cur_loc_id, X0, jumps = pop!(waiting_list)
        cur_loc = HS.modes[cur_loc_id]

        # compute reach tube
        options[:T] = time_horizon - X0.t_start
        reach_tube =
            solve(ContinuousSystem(cur_loc.A, X0.X, cur_loc.U), options)

        # take intersection with source invariant
        reach_tube_in_invariant =
            intersect_reach_tube_invariant(reach_tube.Xk, X0, cur_loc.X, N)
        append!(rset, reach_tube_in_invariant)

        if jumps == max_jumps
            continue
        end

        discrete_post!(waiting_list, HS, cur_loc_id, reach_tube_in_invariant,
            jumps, N)
    end
    options[:T] = time_horizon
    return ReachSolution(rset, options)
end

function intersect_reach_tube_invariant(reach_tube, X0, invariant, N)
    # TODO temporary conversion to HPolytope
    @assert invariant isa HalfSpace
    invariant = HPolytope([invariant])

    # TODO First check for empty intersection, which can be more efficient.
    #      However, we need to make sure that the emptiness check does not just
    #      compute the concrete intersection; otherwise, we would do the work
    #      twice. This is currently the case for 'Polyhedra' polytopes.
    intersections = Vector{ReachSet{LazySet{N}, N}}()
    for reach_set in reach_tube
        rs = reach_set.X
        # TODO temporary workaround for 1D sets
        if dim(rs) == 1
            reach_tube = VPolytope(vertices_list(
                Approximations.overapproximate(rs, LazySets.Interval)))
        # TODO offer a lazy intersection here
        # TODO offer more options instead of taking the VPolytope intersection
        elseif rs isa CartesianProductArray
            reach_tube = VPolytope(vertices_list(rs))
        else
            error("unsupported set type for reach tube: $(typeof(rs))")
        end
        intersect = intersection(invariant, reach_tube)
        if isempty(intersect)
            break
        end
        push!(intersections, ReachSet{LazySet{N}, N}(intersect,
            reach_set.t_start + X0.t_start, reach_set.t_end + X0.t_end))
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
    options.dict[:inout_map] =
        inout_map_reach(options[:partition], options[:blocks], options[:n])

    if options[:project_reachset]
        options[:output_function] = nothing
    else
        options[:output_function] = options[:projection_matrix]
    end

    return options
end
