using HybridSystems
using LazySets.Approximations.box_approximation
export AbstractSolution,
       ReachSolution,
       CheckSolution,
       solve_cont,
       solve_hybrid,
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

"""
    solve(system, options)  or  solve(system, :key1 => val1, [...], keyK => valK)

Solves a reachability problem s.t. the given options.
If some options are not defined, we may fall back to default values.

### Input

- `system`  -- a (discrete or continuoues) system specification
- `options` -- options for solving the problem

### Output

A solution object whose content depends on the input options.

### Notes

To see all available input options, see
`keys(Reachability.available_keywords.dict)`.
"""
function solve_cont(system::AbstractSystem,
               options_input::Options)::AbstractSolution

    # ==========
    # Dimensions
    # ==========
    options_input.dict[:n] = statedim(system)

    # =======
    # Options
    # =======
    # solver-specific options (adds default values for unspecified options)
    options = validate_solver_options_and_add_default_values!(options_input)

    # ===================
    # Time discretization
    # ===================
    if system isa InitialValueProblem{<:AbstractContinuousSystem}
        info("Time discretization...")
        tic()
        Δ = discretize(
            system,
            options[:δ],
            approx_model=options[:approx_model],
            pade_expm=options[:pade_expm],
            lazy_expm=options[:lazy_expm_discretize],
            lazy_sih=options[:lazy_sih]
            )
        tocc()
    else
        Δ = system
    end

    # ==============
    # Transformation
    # ==============
    if options[:coordinate_transformation] != ""
        info("Transformation...")
        tic()
        (Δ, transformation_matrix) = transform(Δ, options[:coordinate_transformation])
        tocc()
    else
        transformation_matrix = nothing
    end

    # ===================================
    # Sparse/dense/lazy matrix conversion
    # ===================================
    A = Δ.s.A
    create_new_system = false
    if !options[:lazy_expm] && options[:lazy_expm_discretize]
        # convert SparseMatrixExp to eplicit matrix
        info("Making lazy matrix exponential explicit...")
        tic()
        n = options.dict[:n]
        if options[:assume_sparse]
            B = sparse(Int[], Int[], eltype(A)[], n, n)
        else
            B = Matrix{eltype(A)}(n, n)
        end
        for i in 1:n
            B[i, :] = get_row(A, i)
        end
        A = B
        create_new_system = true
        tocc()
    end
    if options[:assume_sparse]
        if A isa SparseMatrixExp
            # ignore this case
        elseif !method_exists(sparse, Tuple{typeof(A)})
            info("`assume_sparse` option cannot be applied to a matrix of " *
                 "type $(typeof(A)) and will be ignored")
        elseif !(A isa AbstractSparseMatrix)
            # convert to sparse matrix
            A = sparse(A)
            create_new_system = true
        end
    end
    if create_new_system
        # set new matrix
        if method_exists(inputset, Tuple{typeof(Δ.s)})
            Δ = DiscreteSystem(A, Δ.x0, inputset(Δ))
        else
            Δ = DiscreteSystem(A, Δ.x0)
        end
    end

    if options[:mode] == "reach"
        # ================================
        # Input -> Output variable mapping
        # ================================
        options.dict[:inout_map] =
            inout_map_reach(options[:partition], options[:blocks], options[:n])

        # ============================
        # Reachable states computation
        # ============================
        info("Reachable States Computation...")
        tic()
        output_function = options[:project_reachset] ?
            nothing :
            options[:projection_matrix]
        Rsets = reach(
            Δ,
            options[:N];
            algorithm=options[:algorithm],
            ε_init=options[:ε_init],
            set_type_init=options[:set_type_init],
            ε_iter=options[:ε_iter],
            set_type_iter=options[:set_type_iter],
            assume_sparse=options[:assume_sparse],
            assume_homogeneous=options[:assume_homogeneous],
            lazy_X0=options[:lazy_X0],
            blocks=options[:blocks],
            partition=options[:partition],
            block_types_init=options[:block_types_init],
            block_types_iter=options[:block_types_iter],
            template_directions_init=options[:template_directions_init],
            template_directions_iter=options[:template_directions_iter],
            lazy_inputs_interval=options[:lazy_inputs_interval],
            output_function=output_function
            )
        info("- Total")
        tocc()

        # ==========
        # Projection
        # ==========
        if options[:project_reachset] || options[:projection_matrix] != nothing
            info("Projection...")
            tic()
            RsetsProj = project(Rsets, options;
                                transformation_matrix=transformation_matrix)
            tocc()
            return ReachSolution(RsetsProj, options)
        end

        return ReachSolution(Rsets, options)

    elseif options[:mode] == "check"
        # ============================================
        # Input -> Output variable mapping in property
        # ============================================
        options.dict[:property] = inout_map_property(options[:property],
            options[:partition], options[:blocks], options[:n])

        # =================
        # Property checking
        # =================
        info("Property Checking...")
        tic()
        answer = check_property(
            Δ,
            options[:N];
            algorithm=options[:algorithm],
            ε_init=options[:ε_init],
            set_type_init=options[:set_type_init],
            ε_iter=options[:ε_iter],
            set_type_iter=options[:set_type_iter],
            assume_sparse=options[:assume_sparse],
            assume_homogeneous=options[:assume_homogeneous],
            lazy_X0=options[:lazy_X0],
            blocks=options[:blocks],
            partition=options[:partition],
            block_types_init=options[:block_types_init],
            block_types_iter=options[:block_types_iter],
            template_directions_init=options[:template_directions_init],
            template_directions_iter=options[:template_directions_iter],
            eager_checking=options[:eager_checking],
            property=options[:property],
            lazy_inputs_interval=options[:lazy_inputs_interval]
            )
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
        error("Unsupported mode.")
    end
end

solve_cont(system::AbstractSystem, options::Pair{Symbol,<:Any}...) =
    solve_cont(system, Options(Dict{Symbol,Any}(options)))

"""
    project(Rsets, options; [transformation_matrix])

Projects a sequence of sets according to the settings defined in the options.

### Input

- `Rsets`   -- solution of a reachability problem
- `options` -- options structure
- `transformation_matrix` -- (optional, default: nothing) matrix implementing
                              the transformation)

### Notes

A projection matrix can be given in the options structure, or passed as a
dictionary entry.
"""
function project(Rsets::Vector{<:LazySet}, options::Options;
                 transformation_matrix=nothing)
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
                              transformation_matrix=transformation_matrix,
                              projection_matrix=options[:projection_matrix],
                              output_function=output_function
                             )
end

project(reach_sol::AbstractSolution) = project(reach_sol.Xk, reach_sol.options)

project(Rsets::Vector{<:LazySet}, options::Pair{Symbol,<:Any}...) =
    project(Rsets, Options(Dict{Symbol,Any}(options)))


"""
    solve_hybrid(HS::HybridSystem,
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
function solve_hybrid(HS::HybridSystem,
                      X0::LazySet,
                      options::Options)::AbstractSolution
    @assert isdefined(Main, :Polyhedra) "this algorithm needs the package " *
            "'Polyhedra' to be loaded"

    # set options
    options.dict[:project_reachset] = false

    waiting_list = []
    #TODO get start state. For now we assume that it is the first location
    cur_loc_id = 1
    push!(waiting_list,(cur_loc_id, X0))
    i = 0
    rset = []
    while (!isempty(waiting_list) && i < 15) #TODO add variable for max iteration number
        i += 1
        info("Iteration... $i")
        cur_loc_id, X0 = pop!(waiting_list)
        cur_loc = HS.modes[cur_loc_id]
        S = ContinuousSystem(cur_loc.A, X0, cur_loc.U)

        Rsets = solve_cont(S, options)

        # check intersection with source location
        cur_invariant = cur_loc.X
        intersectedRset =
            intersect_reach_tubes_invariant(Rsets.Xk, cur_invariant)
        push!(rset, intersectedRset)

        # TODO temp assumptions
        # eg. what if we have LazySets.Hyperplane in the array?
        @assert cur_invariant isa HPolytope

        for trans in out_transitions(HS, cur_loc_id)
            info("Considering transition: $trans")
            target_loc_id = target(HS, trans)
            target_loc = HS.modes[target(HS, trans)]
            target_invariant = target_loc.X
            guard = HS.resetmaps[target_loc_id].X
            assignment = HS.resetmaps[target_loc_id].A

            # TODO temp assumptions
            # eg. what if we have LazySets.Hyperplane in the array?
            @assert target_invariant isa HPolytope && guard isa HPolytope

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
                push!(waiting_list,(target(HS, trans), rh))
            end

            # TODO check intersection with forbidden states
            # push!(waiting_list,(target(HS, trans), rsetHull))
        end

        info("End of iteration $i")
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
