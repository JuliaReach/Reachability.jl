__precompile__()
"""
This is the main module and provides interfaces for specifying and solving reachability problems.
"""
module Reachability

include("Systems/Systems.jl")
include("Utils/Utils.jl")
include("ReachSets/ReachSets.jl")
include("Properties/Properties.jl")
include("Transformations/Transformations.jl")

using Reexport
@reexport using LazySets, Reachability.Utils, Reachability.Systems
using Reachability.ReachSets, Reachability.Properties,
      Reachability.Transformations

export Property, Clause

"""
    Options

Type that wraps a dictionary used for options.

FIELDS:

- ``dict`` -- the wrapped dictionary
"""
struct Options
    dict::Dict{Symbol,Any}
    Options(args::Pair{Symbol,<:Any}...) = new(Dict{Symbol,Any}(args))
    Options(dict::Dict{Symbol,Any}) = new(dict)
end


import Base.merge

"""
    merge(op1, opn)

Merges two `Options` objects by just falling back to the wrapped `Dict` fields.
Values are inserted in the order in which the function arguments occur, i.e.,
for conflicting keys a later object overrides a previous value.

INPUT:

- ``op1`` -- first options object
- ``opn`` -- list of options objects
"""
function merge(op1::Options, opn::Options...)::Options
    dict = Dict{Symbol,Any}(op1.dict)
    for i in 1 : length(opn)
        merge!(dict, opn[i].dict)
    end
    return Options(dict)
end


import Base.getindex

"""
    getindex(op, sym)

Returns the value stored for key `sym`.

INPUT:

- ``op`` -- options object
- ``sym`` -- key
"""
function getindex(op::Options, sym::Symbol)
    return getindex(op.dict, sym)
end

export Options, merge, getindex


"""
    validate_solver_options_and_add_default_values!(system, options)

Validates that the given solver options are supported and adds default values for most
unspecified options.

INPUT:

- ``system``  -- a system
- ``options`` -- an `Options` object, a dictionary of options
                 Supported options:
                 - `:mode` (main analysis mode)
                 - `:approx_model` (switch for bloating/continuous time analysis)
                 - `:property` (a safety property)
                 - `:δ` (time step)
                 - `:N` (number of time steps)
                 - `:T` (time horizon)
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
"""
function validate_solver_options_and_add_default_values!(system::System, options::Options)::Options

    dict = options.dict

    # use default values for unspecified options
    # TODO<notification> print message to user for each default option, depending on verbosity level
    if !haskey(dict, :mode)
        dict[:mode] = "reach"
    end
    if !haskey(dict, :approx_model)
        dict[:approx_model] = "forward"
    end
    if !haskey(dict, :property)
        # no default value
    end
    if !haskey(dict, :algorithm)
        dict[:algorithm] = "explicit"
    end
    if !haskey(dict, :blocks)
        dict[:blocks] = [1]
    end
    if !haskey(dict, :iterative_refinement)
        dict[:iterative_refinement] = false
    end
    if !haskey(dict, :ɛ)
        dict[:ɛ] = Inf
    end
    if !haskey(dict, :lazy_expm)
        dict[:lazy_expm] = false
    end
    if !haskey(dict, :assume_sparse)
        dict[:assume_sparse] = false
    end
    if !haskey(dict, :pade_expm)
        dict[:pade_expm] = false
    end
    if !haskey(dict, :lazy_X0)
        dict[:lazy_X0] = false
    end
    if !haskey(dict, :discr_algorithm)
        dict[:discr_algorithm] = "forward"
    end
    if !haskey(dict, :coordinate_transformation)
        dict[:coordinate_transformation] = ""
    end
    if !haskey(dict, :assume_homogeneous)
        dict[:assume_homogeneous] = false
    end
    if !haskey(dict, :projection_matrix)
        dict[:projection_matrix] = nothing
    end
    if !haskey(dict, :apply_projection)
        dict[:apply_projection] = true
    end
    if !haskey(dict, :plot_vars)
        dict[:plot_vars] = [0, 1]
    end
    
    # special options: δ, N, T
    check_and_add_δ_N_T!(dict)

    # validation
    for kv_pair in dict
        key::Symbol = kv_pair.first
        value = kv_pair.second

        # define type/domain constraints for each known key
        domain_constraints = (v  ->  true)
        if key == :mode
            expected_type = String
            domain_constraints = (v::String  ->  v in ["reach", "check"])
            if value == "check" && !haskey(dict, :property)
                error("No property has been defined.")
            end
        elseif key == :approx_model
            expected_type = String
        elseif key == :property
            expected_type = Property
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

    return options
end


"""
Handling of the special trio `:δ`, `:N`, `:T`.
Usually two of them should be defined and the third one is automatically inferred.
If three of them are defined, they must be consistent.
If only `:T` is defined, we use `:N = 100`.

INPUT:

- ``dict`` -- a dictionary of options
"""
function check_and_add_δ_N_T!(dict::Dict{Symbol,Any})
    defined = 0
    if haskey(dict, :δ)
        value = dict[:δ]
        if !(value isa Float64 && value > 0.)
            error(string(value) * " is not a valid value for option :δ.")
        end
        defined += 1
    end
    if haskey(dict, :N)
        value = dict[:N]
        if !(value isa Int64 && value > 0)
            error(string(value) * " is not a valid value for option :N.")
        end
        defined += 1
    end
    if haskey(dict, :T)
        value = dict[:T]
        if !(value isa Float64 && value > 0.)
            error(string(value) * " is not a valid value for option :T.")
        end
        defined += 1
    end

    if defined == 1 && haskey(dict, :T)
        dict[:N] = 100
        dict[:δ] = dict[:T] / dict[:N]
    elseif defined == 2
        if !haskey(dict, :δ)
            dict[:δ] = dict[:T] / dict[:N]
        end
        if !haskey(dict, :N)
            dict[:N] = (Int64)(ceil(dict[:T] / dict[:δ]))
        end
        if !haskey(dict, :T)
            dict[:T] = dict[:N] * dict[:δ]
        end
    elseif defined == 3
        N_computed = (Int64)(ceil(dict[:T] / dict[:δ]))
        if N_computed != dict[:N]
            error("Values for :δ, :N, and :T were defined, but they are inconsistent.")
        end
    else
        error("Need two of the following options: :δ, :N, :T")
    end
end


"""
    solve(system, options)  or  solve(system, :key1 => val1, [...], keyK => valK)

Solves a reachability problem s.t. the given options.
If some options are not defined, we may fall back to default values.

INPUT:

- ``system`` -- a (discrete or continuoues) system specification
- ``options`` -- options for solving the problem
"""
function solve(system::Union{ContinuousSystem, DiscreteSystem}, options::Options)

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
    # solver-specific options (uses default for unspecified options)
    validate_solver_options_and_add_default_values!(system, options)

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
            return RsetsProj
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

export solve


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


project(Rsets::Union{Vector{CartesianProductArray}, Vector{HPolygon}}, n::Int64,
        options::Pair{Symbol,<:Any}...) =
            project(Rsets, n, Options(Dict{Symbol,Any}(options)))

export project


"""
    validate_plot_options_and_add_default_values!(polygons, options)

Validates that the given plot options are supported and adds default values for all
unspecified options.

INPUT:

- ``polygons`` -- sequence of polygons
- ``options``  -- an `Options` object, a dictionary of options
                  Supported options:
                  - `:plot_indices` (which indices of the sequence to plot)
                  - `:plot_backend` (which backend should be used)
                  - `:gridlines` (should grid lines be used)
                  - `:plot_name` (name of the plot file)
                  - `:plot_vars` (variables to plot)
                  - `:plot_labels` (axis labels)
                  - `:plot_color` (plot color)
"""
function validate_plot_options_and_add_default_values!(polygons::Vector{HPolygon}, options::Options)
    dict = options.dict

    # use default values for unspecified options
    # TODO<notification> print message to user for each default option, depending on verbosity level
    if !haskey(dict, :plot_indices)
        dict[:plot_indices] = 1:length(polygons)
    end
    if !haskey(dict, :plot_backend)
        dict[:plot_backend] = "pyplot_savefig"
    end
    if !haskey(dict, :gridlines)
        dict[:gridlines] = true
    end
    if !haskey(dict, :plot_name)
        dict[:plot_name] = "plot.png"
    end
    if !haskey(dict, :plot_vars)
        dict[:plot_vars] = [0, 1]
    end
    if !haskey(dict, :plot_labels)
        dict[:plot_labels] = add_plot_labels(dict[:plot_vars])
    end
    if !haskey(dict, :plot_color)
        dict[:plot_color] = "blue"
    end

    # validation
    for kv_pair in dict
        key::Symbol = kv_pair.first

        # define type/domain constraints for each known key
	domain_constraints = (v  ->  true)
	if key == :plot_indices
            expected_type = AbstractVector{Int64}
        elseif key == :plot_backend
            expected_type = String
            domain_constraints = (v::String  ->  v in ["", "pyplot_savefig", "pyplot_inline", "gadfly"])
        elseif key == :gridlines
            expected_type = Bool
        elseif key == :plot_name
            expected_type = String
        elseif key == :plot_labels
            expected_type = Vector{String}
            domain_constraints = (v::Vector{String}  ->  length(v) == 2)
        elseif key == :plot_vars
            expected_type = Vector{Int64}
            domain_constraints = (v::Vector{Int64}  ->  length(v) == 2)
        elseif key == :plot_color
            expected_type = String
        else
            error("Unrecognized option '" * string(key) * "' found.")
        end

        value = kv_pair.second
        # check value type
        if !(value isa expected_type)
            error("Option :" * string(key) * " must be of '" * string(expected_type) * "' type.")
        end
        # check value domain constraints
        if !domain_constraints(value)
            error(string(value) * " is not a valid value for option " * string(key) * ".")
        end
    end

    return options
end


"""
    plot(polygons, [options])

Plots a sequence of polygons with the given options.

INPUT:

- ``polygons`` -- sequence of polygons
- ``options`` -- options for plotting
"""
function plot(polygons::Vector{HPolygon}, options::Options)
    options = validate_plot_options_and_add_default_values!(polygons, options)

    plot_polygon(
        view(polygons, options[:plot_indices]),
        backend=options[:plot_backend],
        gridlines=options[:gridlines],
        name=options[:plot_name],
        plot_labels=options[:plot_labels],
        color=options[:plot_color]
        )
    nothing
end


plot(polygons::Vector{HPolygon}, options::Pair{Symbol,<:Any}...) = plot(polygons, Options(Dict{Symbol,Any}(options)))

export plot

end # module
