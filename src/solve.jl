export solve

function default_operator(system::InitialValueProblem)
    s = system.s
    if s isa LinearContinuousSystem ||
            s isa LinearControlContinuousSystem ||
            s isa ConstrainedLinearContinuousSystem ||
            s isa ConstrainedLinearControlContinuousSystem ||
            s isa LinearDiscreteSystem ||
            s isa LinearControlDiscreteSystem ||
            s isa ConstrainedLinearDiscreteSystem ||
            s isa ConstrainedLinearControlDiscreteSystem
        op = BFFPSV18()
    else
        error("no default reachability algorithm available for system of " *
              "type $(typeof(system))")
    end
    return op
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
               op::PostOperator=default_operator(system))
    solve!(system, Options(copy(options.dict)), op=op)
end

solve(system::AbstractSystem, options::Pair{Symbol,<:Any}...) =
    solve(system, Options(Dict{Symbol,Any}(options)))

function solve!(system::InitialValueProblem{<:SYS},
                options_input::Options;
                op::PostOperator=default_operator(system)
               )::AbstractSolution where {SYS<:Union{AbstractContinuousSystem,
                                                     AbstractDiscreteSystem}}
    options = init(op, system, options_input)

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

    post(op, system, options)
end

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
               options_input::Options)::AbstractSolution where {N}
    @assert isdefined(Main, :Polyhedra) "this algorithm needs the package " *
            "'Polyhedra' to be loaded"

    options, time_horizon, max_jumps, delete_N = init_hybrid!(HS, options_input)

    #TODO get initial location. For now we assume that it is the first location
    # waiting_list entries:
    # - (discrete) location
    # - (set of) continuous-time reach sets
    # - number of previous jumps
    waiting_list = [(1, ReachSet{LazySet{N}, N}(X0, zero(N), zero(N)), 0)]
    # passed_list maps the (discrete) location to the (set of) continuous-time
    # reach sets
    passed_list = Vector{Vector{ReachSet{LazySet{N}, N}}}(nstates(HS))
    Rsets = Vector{ReachSet{LazySet{N}, N}}()
    while (!isempty(waiting_list))
        cur_loc_id, X0, jumps = pop!(waiting_list)
        cur_loc = HS.modes[cur_loc_id]

        # compute reach tube
        options_copy = Options(copy(options.dict))
        options_copy.dict[:T] = time_horizon - X0.t_start
        options_copy.dict[:project_reachset] = false
        delete!(options_copy.dict, :inout_map)
        if delete_N # TODO add more conditions or fix option clashes in general
            delete!(options_copy.dict, :N)
        end
        if haskey(options_copy.dict, :block_types) &&
                options_copy.dict[:block_types] == nothing
            delete!(options_copy.dict, :block_types)
        end
        if haskey(options_copy.dict, :blocks)
            delete!(options_copy.dict, :blocks)
        end
        reach_tube =
            solve(ContinuousSystem(cur_loc.A, X0.X, cur_loc.U), options_copy)

        # take intersection with source invariant
        reach_tube_in_invariant =
            intersect_reach_tube_invariant(reach_tube.Xk, X0, cur_loc.X, N)
        append!(Rsets, reach_tube_in_invariant)

        if jumps == max_jumps
            continue
        end

        discrete_post!(waiting_list, passed_list, HS, cur_loc_id,
            reach_tube_in_invariant, jumps, N)
    end

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

function init_hybrid!(system, options_input)
    delete_N = !haskey(options_input.dict, :N)
    options_input.dict[:n] = statedim(system, 1)

    # solver-specific options (adds default values for unspecified options)
    options = validate_solver_options_and_add_default_values!(options_input)

    # Input -> Output variable mapping
    options.dict[:inout_map] =
        inout_map_reach(options[:partition], options[:blocks], options[:n])

    time_horizon = options[:T]
    max_jumps = options[:max_jumps]

    return options, time_horizon, max_jumps, delete_N
end
