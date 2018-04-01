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
function solve(system::AbstractSystem,
               options_input::Options)::AbstractSolution

    # ==========
    # Dimensions
    # ==========
    dimension = statedim(system)
    if isodd(dimension)
         system = add_dimension(system)
    end
    options_input.dict[:n] = dimension
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
            lazy_expm=options[:lazy_expm],
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

    # ==============================
    # Sparse/dense matrix conversion
    # ==============================
    A = Δ.s.A
    if options[:assume_sparse]
        if A isa SparseMatrixExp || !method_exists(sparse, Tuple{typeof(Δ.s.A)})
            info("`assume_sparse` option cannot be applied to a matrix of type $(typeof(Δ.s.A)) and will be ignored")
        else
            A = sparse(Δ.s.A)
        end
        Δ = DiscreteSystem(A, Δ.x0, inputset(Δ))
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
            block_types_iter=options[:block_types_iter]
            )
        info("- Total")
        tocc()

        # ==========
        # Projection
        # ==========
        if options[:apply_projection]
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
            property=options[:property]
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

solve(system::AbstractSystem, options::Pair{Symbol,<:Any}...) =
    solve(system, Options(Dict{Symbol,Any}(options)))

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
    RsetsProj = project_reach(plot_vars,
                              reduced_n,
                              options[:δ],
                              Rsets,
                              options[:algorithm],
                              ε=options[:ε_proj],
                              set_type=options[:set_type_proj],
                              transformation_matrix=transformation_matrix,
                              projection_matrix=options[:projection_matrix]
                             )
end

project(reach_sol::AbstractSolution) = project(reach_sol.Xk, reach_sol.options)

project(Rsets::Vector{<:LazySet}, options::Pair{Symbol,<:Any}...) =
    project(Rsets, Options(Dict{Symbol,Any}(options)))
