export solve, project

"""
    validate_solver_options_and_add_default_values!(options)

Validates that the given solver options are supported and adds default values for most
unspecified options.

INPUT:

- ``options`` -- an `Options` object, a dictionary of options

Supported options:
- `:mode` (main analysis mode)
- `:approx_model` (switch for bloating/continuous time analysis)
- `:property` (a safety property)
- `:δ` (time step)
  - alias: `:sampling_time`
- `:N` (number of time steps)
- `:T` (time horizon)
  - alias: `:time_horizon`
- `:algorithm` (algorithm backend)
- `:blocks` (blocks of interest)
- `:iterative_refinement` (switch for refining precision/directions)
- `:ɛ` (precision threshold, s. :iterative_refinement)
- `:lazy_expm` (lazy matrix exponential)
- `:assume_sparse` (switch for sparse matrices)
- `:pade_expm` (switch for using Pade approximant method)
- `:lazy_X0` (switch for keeping the initial states a lazy set)
- `:discr_algorithm` (discretization algorithm)
- `:coordinate_transformation` (coordinate transformation method)
- `:assume_homogeneous` (switch for ignoring inputs)
- `:projection_matrix` (projection matrix)
- `:apply_projection` (switch for applying projection)
- `:plot_vars` (variables for projection and plotting)
  - alias: `:output_variables`

We add default values for almost all undefined options, i.e., modify the input
options. The benefit is that the user can print the options that were actually
used. For supporting aliases, we create another copy that is actually used where
we only keep those option names that are used internally.
"""
function validate_solver_options_and_add_default_values!(options::Options)::Options
    dict = options.dict
    options_copy = Options()
    dict_copy = options_copy.dict

    # check for aliases and use default values for unspecified options
    check_aliases_and_add_default_value!(dict, dict_copy, [:mode], "reach")
    check_aliases_and_add_default_value!(dict, dict_copy, [:approx_model], "forward")
    check_aliases_and_add_default_value!(dict, dict_copy, [:property], nothing)
    check_aliases_and_add_default_value!(dict, dict_copy, [:algorithm], "explicit")
    check_aliases_and_add_default_value!(dict, dict_copy, [:blocks], [1])
    check_aliases_and_add_default_value!(dict, dict_copy, [:iterative_refinement], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:ɛ], Inf)
    check_aliases_and_add_default_value!(dict, dict_copy, [:lazy_expm], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:assume_sparse], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:pade_expm], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:lazy_X0], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:discr_algorithm], "forward")
    check_aliases_and_add_default_value!(dict, dict_copy, [:coordinate_transformation], "")
    check_aliases_and_add_default_value!(dict, dict_copy, [:assume_homogeneous], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:projection_matrix], nothing)
    check_aliases_and_add_default_value!(dict, dict_copy, [:apply_projection], true)
    check_aliases_and_add_default_value!(dict, dict_copy, [:plot_vars, :output_variables], [0, 1])

    # special options: δ, N, T
    check_and_add_δ_N_T!(dict, dict_copy)

    # validation
    for kv_pair in dict_copy
        key::Symbol = kv_pair.first
        value = kv_pair.second

        # define type/domain constraints for each known key
        domain_constraints = (v  ->  true)
        if key == :mode
            expected_type = String
            domain_constraints = (v::String  ->  v in ["reach", "check"])
            if value == "check" && dict_copy[:property] == nothing
                error("No property has been defined.")
            end
        elseif key == :approx_model
            expected_type = String
        elseif key == :property
            expected_type = Union{Property, Void}
        elseif key == :δ
            expected_type = Float64
            domain_constraints = (v::Float64  ->  v > 0.)
        elseif key == :N
            expected_type = Int64
            domain_constraints = (v::Int64  ->  v > 0)
        elseif key == :T
            expected_type = Float64
            domain_constraints = (v::Float64  ->  v > 0.)
        elseif key == :algorithm
            expected_type = String
            domain_constraints = (v::String  ->  v in ["explicit"])
        elseif key == :blocks
            expected_type = AbstractVector{Int64}
        elseif key == :iterative_refinement
            expected_type = Bool
        elseif key == :ɛ
            expected_type = Float64
            domain_constraints = (v::Float64  ->  v > 0.)
        elseif key == :lazy_expm
            expected_type = Bool
        elseif key == :assume_sparse
            expected_type = Bool
        elseif key == :pade_expm
            expected_type = Bool
        elseif key == :lazy_X0
            expected_type = Bool
        elseif key == :discr_algorithm
            expected_type = String
            domain_constraints = (v::String  ->  v in ["forward"])
        elseif key == :coordinate_transformation
            expected_type = String
            domain_constraints = (v::String  ->  v in ["", "schur"])
        elseif key == :assume_homogeneous
            expected_type = Bool
        elseif key == :projection_matrix
            expected_type = Union{SparseMatrixCSC{Float64, Int64}, Void}
        elseif key == :apply_projection
            expected_type = Bool
        elseif key == :plot_vars
            expected_type = Vector{Int64}
            domain_constraints = (v::Vector{Int64}  ->  length(v) == 2)
        else
            error("Unrecognized option '" * string(key) * "' found.")
        end

        # check value type
        if !(value isa expected_type)
            error("Option :" * string(key) * " must be of '" * string(expected_type) * "' type.")
        end
        # check value domain constraints
        if !domain_constraints(value)
            error(string(value) * " is not a valid value for option " * string(key) * ".")
        end
    end

    return options_copy
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
