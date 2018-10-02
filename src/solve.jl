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
function solve(system::InitialValueProblem, options::Options; algorithm::String=default_algorithm(system))
    solve!(system, Options(copy(options.dict)), algorithm=algorithm)
end

solve(system::AbstractSystem, options::Pair{Symbol,<:Any}...) =
    solve(system, Options(Dict{Symbol,Any}(options)))

function solve!(system::InitialValueProblem, options::Options;
                algorithm::String=default_algorithm(system))::AbstractSolution
    if algorithm == "BFFPSV18"
        options, output_function = init_BFFPSV18!(system, options)

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
            Rsets = reach(
                system,
                options[:N],
                options;
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
            answer = check_property(
                system,
                options[:N],
                options;
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
        output_function = nothing
    else
        output_function = options[:projection_matrix]
    end

    return options, output_function
end

function matrix_conversion(Δ, options)

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
    return Δ
end
