export Rset2DSolution, solve, project

"""
    Rset2DSolution

Type that wraps the two-dimensional reach sets of the solution of a
reachability problem.

FIELDS:

- ``polygons`` -- the list of reachable states, given as polygons in constraint
                  representation
- ``options``  -- the dictionary of options
"""
struct Rsets2DSolution
  polygons::Vector{HPolygon}
  options::Options

  Rsets2DSolution(polygons, args::Options) = new(polygons, args)
  Rset2sDSolution(polygons) = new(polygons, Options())
end

"""
    solve(system, options)  or  solve(system, :key1 => val1, [...], keyK => valK)

Solves a reachability problem s.t. the given options.
If some options are not defined, we may fall back to default values.

INPUT:

- ``system`` -- a (discrete or continuoues) system specification
- ``options`` -- options for solving the problem
"""
function solve(system::Union{ContinuousSystem, DiscreteSystem}, options_input::Options)

    # ==========
    # Dimensions
    # ==========
    dimension = Systems.dim(system)
    if isodd(dimension)
         system = add_dimension(system)
    end

    # =======
    # Options
    # =======
    # solver-specific options (adds default values for unspecified options)
    options = validate_solver_options_and_add_default_values!(options_input)

    # ===================
    # Time discretization
    # ===================
    if system isa ContinuousSystem
        println("Time discretization...")
        tic()
        Δ = discretize(
            system,
            options[:δ],
            approx_model=options[:approx_model],
            pade_expm=options[:pade_expm],
            lazy_expm=options[:lazy_expm]
            )
        toc()
    else
        Δ = system
    end

    # ==============
    # Transformation
    # ==============
    if options[:coordinate_transformation] != ""
        println("Transformation...")
        tic()
        (Δ, transformation_matrix) = transform(
            options[:coordinate_transformation],
            Δ,
            options[:plot_vars]
            )
        toc()
    else
        transformation_matrix = nothing
    end

    if options[:mode] == "reach"
        # ============================
        # Reachable states computation
        # ============================
        println("Reachable States Computation...")
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
        toc()

        # ==========
        # Projection
        # ==========
        if options[:apply_projection]
            println("Projection...")
            tic()
            RsetsProj = project(Rsets, size(Δ.A, 1), options, transformation_matrix)
            toc()
            return Rsets2DSolution(RsetsProj, options)
        end

        return Rsets

    elseif options[:mode] == "check"
        # =================
        # Property checking
        # =================
        println("Property Checking...")
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
        toc()

        if answer == 0
            println("The property is satisfied!")
            return true
        else
            println("The property may be violated at index ", answer, " (time point ", answer * options[:δ], ")!")
            return false
        end
    else
        error("Unsupported mode.")
    end
end

solve(system::Union{ContinuousSystem, DiscreteSystem}, options::Pair{Symbol,<:Any}...) = solve(system, Options(Dict{Symbol,Any}(options)))

"""
    project(Rsets, n, options, [transformation_matrix])

Projects a sequence of sets according to the settings defined in the options.

INPUT:

- ``Rsets``                 -- sequence of sets
- ``n``                     -- dimension
- ``options``               -- options structure
- ``transformation_matrix`` -- (optional, default: `nothing`) transformation matrix
"""
function project(Rsets::Union{Vector{CartesianProductArray}, Vector{HPolygon}},
                 n::Int64, options::Options, transformation_matrix=nothing)

    RsetsProj = project_reach(
        options[:plot_vars],
        n,
        options[:δ],
        Rsets,
        options[:algorithm],
        ɛ=options[:ɛ],
        transformation_matrix=transformation_matrix,
        projection_matrix=options[:projection_matrix]
        )
end

project(Rsets::Union{Vector{CartesianProductArray}, Vector{HPolygon}},
        n::Int64,
        options::Pair{Symbol,<:Any}...) = project(Rsets, n, Options(Dict{Symbol,Any}(options)))
