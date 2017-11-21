export ReachSolution, AbstractReachSolution, solve, project

"""
    AbstractReachSolution

Abstract type representing the solution of a rechability problem.
"""
abstract type AbstractReachSolution end

"""
    ReachSolution

Type that wraps a the solution of a reachability problem as a sequence of
lazy sets, and a dictionary of options.

### Fields

- `Xk`       -- the list of reachable states
- `options`  -- the dictionary of options

### Note

If the solution has been projected in 2D, the sequence `Xk` is an array
of polygons in constraint representation. In high-dimensions this is a sequence
of cartesian product arrays of low-dimensional sets.
"""
struct ReachSolution{S<:LazySet} <: AbstractReachSolution
  Xk::Vector{S}
  options::Options
end
ReachSolution(Xk::Vector{S}, options) where {S<:LazySet} = ReachSolution{S}(Xk, options)
ReachSolution(Xk::Vector{S}) where {S<:LazySet} = ReachSolution{S}(Xk, Options())

"""
    solve(system, options)  or  solve(system, :key1 => val1, [...], keyK => valK)

Solves a reachability problem s.t. the given options.
If some options are not defined, we may fall back to default values.

### Input

- `system`  -- a (discrete or continuoues) system specification
- `options` -- options for solving the problem
"""
function solve(system::AbstractSystem,
               options_input::Options)::AbstractReachSolution

    # ==========
    # Dimensions
    # ==========
    dimension = Systems.dim(system)
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
    if system isa ContinuousSystem
        info("Time discretization...")
        tic()
        Δ = discretize(
            system,
            options[:δ],
            approx_model=options[:approx_model],
            pade_expm=options[:pade_expm],
            lazy_expm=options[:lazy_expm]
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
        (Δ, transformation_matrix) = transform(
            options[:coordinate_transformation],
            Δ,
            options[:plot_vars]
            )
        tocc()
    else
        transformation_matrix = nothing
    end

    if options[:mode] == "reach"
        # ============================
        # Reachable states computation
        # ============================
        info("Reachable States Computation...")
        tic()
        Rsets = reach(
            Δ,
            options[:N];
            algorithm=options[:algorithm],
            ɛ=options[:ɛ],
            blocks=options[:blocks],
            assume_sparse=options[:assume_sparse],
            iterative_refinement=options[:iterative_refinement],
            assume_homogeneous=options[:assume_homogeneous],
            lazy_X0=options[:lazy_X0]
            )
        tocc()

        # ==========
        # Projection
        # ==========
        if options[:apply_projection]
            info("Projection...")
            tic()
            RsetsProj = project_reach(options[:plot_vars],
                                      dimension,
                                      options[:δ],
                                      Rsets,
                                      options[:algorithm];
                                      ɛ=options[:ɛ],
                                      transformation_matrix=transformation_matrix,
                                      projection_matrix=options[:projection_matrix])
            tocc()
            return ReachSolution(RsetsProj, options)
        end

        return ReachSolution(Rsets, options)

    elseif options[:mode] == "check"
        # =================
        # Property checking
        # =================
        info("Property Checking...")
        tic()
        answer = check_property(
            Δ,
            options[:N];
            algorithm=options[:algorithm],
            ɛ=options[:ɛ],
            blocks=options[:blocks],
            assume_sparse=options[:assume_sparse],
            iterative_refinement=options[:iterative_refinement],
            assume_homogeneous=options[:assume_homogeneous],
            property=options[:property]
            )
        tocc()

        if answer == 0
            info("The property is satisfied!")
            return true
        else
            info("The property may be violated at index $answer, (time point $(answer * options[:δ]))!")
            return false
        end
    else
        error("Unsupported mode.")
    end
end

solve(system::Union{ContinuousSystem, DiscreteSystem}, options::Pair{Symbol,<:Any}...) = solve(system, Options(Dict{Symbol,Any}(options)))

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
function project(Rsets::Union{Vector{CartesianProductArray}, Vector{HPolygon}},
                 options::Options; transformation_matrix=nothing)

    RsetsProj = project_reach(options[:plot_vars],
                              options[:n],
                              options[:δ],
                              Rsets,
                              options[:algorithm],
                              ɛ=options[:ɛ],
                              transformation_matrix=transformation_matrix,
                              projection_matrix=options[:projection_matrix]
                              )
end

project(reach_sol::AbstractReachSolution) = project(reach_sol.Xk, reach_sol.options)

project(Rsets::Union{Vector{CartesianProductArray}, Vector{HPolygon}},
        options::Pair{Symbol,<:Any}...) = project(Rsets, Options(Dict{Symbol,Any}(options)))
