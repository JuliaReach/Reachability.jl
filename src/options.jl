import Base: merge, getindex

export Options, merge, getindex

available_keywords = Set{Symbol}([])

"""
    Options

Type that wraps a dictionary used for options.

FIELDS:

- `dict` -- the wrapped dictionary
"""
struct Options
    dict::Dict{Symbol,Any}
    Options(args::Pair{Symbol,<:Any}...) = new(Dict{Symbol,Any}(args))
    Options(dict::Dict{Symbol,Any}) = new(dict)
end

"""
    merge(op1, opn)

Merges two `Options` objects by just falling back to the wrapped `Dict` fields.
Values are inserted in the order in which the function arguments occur, i.e.,
for conflicting keys a later object overrides a previous value.

INPUT:

- `op1` -- first options object
- `opn` -- list of options objects
"""
function merge(op1::Options, opn::Options...)::Options
    dict = Dict{Symbol,Any}(op1.dict)
    for i in 1 : length(opn)
        merge!(dict, opn[i].dict)
    end
    return Options(dict)
end

"""
    getindex(op, sym)

Returns the value stored for key `sym`.

INPUT:

- `op` -- options object
- `sym` -- key
"""
function getindex(op::Options, sym::Symbol)
    return getindex(op.dict, sym)
end


"""
    validate_solver_options_and_add_default_values!(options)

Validates that the given solver options are supported and adds default values for most
unspecified options.

### Input

- `options` -- an `Options` object, a dictionary of options

Supported options:

- `:verbosity`     -- controls logging output
- `:mode`          -- main analysis mode
- `:approx_model`  -- switch for bloating/continuous time analysis
- `:property`      -- a safety property
- `:δ`             -- time step; alias: `:sampling_time`
- `:N`             -- number of time steps
- `:T`             -- time horizon; alias `:time_horizon`
- `:algorithm`     -- algorithm backend
- `:n`             -- system's dimension
- `:blocks`        -- blocks of interest
- `:iterative_refinement` -- switch for refining precision/directions
- `:ɛ`             -- precision threshold, see also: `:iterative_refinement`
- `:lazy_expm`     -- lazy matrix exponential
- `:assume_sparse` -- switch for sparse matrices
- `:pade_expm`     -- switch for using Pade approximant method
- `:set_type`      -- set type for overapproximation
- `:lazy_X0`       -- switch for keeping the initial states a lazy set
- `:approx_model`  -- approximation model
- `:coordinate_transformation` -- coordinate transformation method
- `:assume_homogeneous`        -- switch for ignoring inputs
- `:projection_matrix`         -- projection matrix
- `:apply_projection`          -- switch for applying projection
- `:plot_vars`                 -- variables for projection and plotting;
                                  alias: `:output_variables`

We add default values for almost all undefined options, i.e., modify the input
options. The benefit is that the user can print the options that were actually
used. For supporting aliases, we create another copy that is actually used where
we only keep those option names that are used internally.
"""
function validate_solver_options_and_add_default_values!(options::Options)::Options
    global LOGGER

    dict = options.dict
    options_copy = Options()
    dict_copy = options_copy.dict

    # first read the verbosity option and set global log level accordingly
    if haskey(dict, :verbosity)
        LOGGER = configure_logger(dict[:verbosity])
        dict_copy[:verbosity] = dict[:verbosity]
    else
        LOGGER = configure_logger()
    end

    # check for aliases and use default values for unspecified options
    check_aliases_and_add_default_value!(dict, dict_copy, [:mode], "reach")
    check_aliases_and_add_default_value!(dict, dict_copy, [:approx_model], "forward")
    check_aliases_and_add_default_value!(dict, dict_copy, [:property], nothing)
    check_aliases_and_add_default_value!(dict, dict_copy, [:algorithm], "explicit")
    check_aliases_and_add_default_value!(dict, dict_copy, [:n], nothing)
    check_aliases_and_add_default_value!(dict, dict_copy, [:blocks], [1])
    check_aliases_and_add_default_value!(dict, dict_copy, [:iterative_refinement], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:ɛ], Inf)
    check_aliases_and_add_default_value!(dict, dict_copy, [:lazy_expm], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:assume_sparse], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:pade_expm], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:set_type], HPolygon)
    check_aliases_and_add_default_value!(dict, dict_copy, [:lazy_X0], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:coordinate_transformation], "")
    check_aliases_and_add_default_value!(dict, dict_copy, [:assume_homogeneous], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:projection_matrix], nothing)
    check_aliases_and_add_default_value!(dict, dict_copy, [:apply_projection], true)
    check_aliases_and_add_default_value!(dict, dict_copy, [:plot_vars, :output_variables], [0, 1])

    # special options: δ, N, T
    check_and_add_δ_N_T!(dict, dict_copy)

    # validate that all input keywords are recognized
    check_valid_option_keywords(dict)

    # validation
    for kv_pair in dict_copy
        key::Symbol = kv_pair.first
        value = kv_pair.second

        # define type/domain constraints for each known key
        domain_constraints = (v  ->  true)
        if key == :verbosity
            expected_type = Union{String, Int}
        elseif key == :mode
            expected_type = String
            domain_constraints = (v::String  ->  v in ["reach", "check"])
            if value == "check" && dict_copy[:property] == nothing
                error("No property has been defined.")
            end
        elseif key == :approx_model
            expected_type = String
            domain_constraints = (v::String  ->  v in ["forward", "backward",
                                                       "firstorder", "nobloating"])
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
        elseif key == :n
            expected_type = Int
            domain_constraints = (v::Int  ->  v > 0)
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
        elseif key == :set_type
            expected_type = Union{Type{HPolygon}, Type{Hyperrectangle}}
        elseif key == :lazy_X0
            expected_type = Bool
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
            error(get_unrecognized_key_message(key))
        end

        # check value type
        if !(value isa expected_type)
            error("Option :$key must be of '$expected_type' type.")
        end
        # check value domain constraints
        if !domain_constraints(value)
            error("$value is not a valid value for option :$key.")
        end
    end

    return options_copy
end

"""
Handling of the special trio `:δ`, `:N`, `:T`.
Usually two of them should be defined and the third one is automatically inferred.
If three of them are defined, they must be consistent.
If only `:T` is defined, we use `:N = 100`.

INPUT:

- `dict`      -- a dictionary of options
- `dict_copy` -- a copy of the dictionary of options for internal names

Supported options:
- `:δ` (time step)
- alias: `:sampling_time`
- `:N` (number of time steps)
- `:T` (time horizon)
- alias: `:time_horizon`
"""
function check_and_add_δ_N_T!(dict::Dict{Symbol,Any}, dict_copy::Dict{Symbol,Any})
    δ_aliases = [:δ, :sampling_time]
    T_aliases = [:T, :time_horizon]
    check_aliases!(dict, dict_copy, δ_aliases)
    check_aliases!(dict, dict_copy, [:N])
    check_aliases!(dict, dict_copy, T_aliases)

    defined = 0
    if haskey(dict_copy, :δ)
        value = dict_copy[:δ]
        if !(value isa Float64 && value > 0.)
            error("$value is not a valid value for option $δ_aliases.")
        end
        defined += 1
    end
    if haskey(dict_copy, :N)
        value = dict_copy[:N]
        if !(value isa Int64 && value > 0)
            error("$value is not a valid value for option :N.")
        end
        defined += 1
    end
    if haskey(dict_copy, :T)
        value = dict_copy[:T]
        if !(value isa Float64 && value > 0.)
            error("$value is not a valid value for option $T_aliases.")
        end
        defined += 1
    end

    if defined == 1 && haskey(dict_copy, :T)
        N = 100
        dict[:N] = N
        dict_copy[:N] = N
        δ = dict_copy[:T] / dict_copy[:N]
        dict[:δ] = δ
        dict_copy[:δ] = δ
    elseif defined == 2
        if !haskey(dict_copy, :δ)
            δ = dict_copy[:T] / dict_copy[:N]
            dict[:δ] = δ
            dict_copy[:δ] = δ
        end
        if !haskey(dict_copy, :N)
            N = (Int64)(ceil(dict_copy[:T] / dict_copy[:δ]))
            dict[:N] = N
            dict_copy[:N] = N
        end
        if !haskey(dict_copy, :T)
            T = dict_copy[:N] * dict_copy[:δ]
            dict[:T] = T
            dict_copy[:T] = T
        end
    elseif defined == 3
        N_computed = (Int64)(ceil(dict_copy[:T] / dict_copy[:δ]))
        if N_computed != dict_copy[:N]
            error("Values for :δ ($(dict_copy[:δ])), :N ($(dict_copy[:N])), and :T ($(dict_copy[:T])) were defined, but they are inconsistent.")
        end
    else
        error("Need two of the following options: :δ, :N, :T (or at least :T and using default values)")
    end
end


"""
    check_aliases!(dict, dict_copy, aliases)

This function has several purposes:

- translate aliases to the option that is used internally
- check that not several aliases were used at the same time

INPUT:

- `dict` -- a dictionary of options
- `dict_copy` -- a copy of the dictionary of options for internal names
- `aliases` -- option aliases; the first name is the one we use internally
"""
function check_aliases!(dict::Dict{Symbol,Any}, dict_copy::Dict{Symbol,Any}, aliases::Vector{Symbol})
    # find aliases and check consistency in case several aliases have been used
    print_warning::Bool = false
    for alias in aliases
        # add alias to global keyword set
        push!(available_keywords, alias)

        if haskey(dict, alias)
            # alias was used
            if haskey(dict_copy, aliases[1])
                # several aliases were used
                if (dict_copy[aliases[1]] != dict[alias])
                    error("Several option aliases were used with different values for aliases $aliases.")
                else
                    print_warning = true
                end
            else
                dict_copy[aliases[1]] = dict[alias]
            end
        end
    end
    if print_warning
        warn("Several option aliases were used for aliases $aliases.")
    end
end

"""
    check_aliases_and_add_default_value!(dict, dict_copy, aliases, default_value, [modify_dict])

This function has several purposes:
- translate aliases to the option that is used internally
- check that not several aliases were used at the same time
- assign default values for undefined options

INPUT:

- `dict` -- a dictionary of options
- `dict_copy` -- a copy of the dictionary of options for internal names
- `aliases` -- option aliases; the first name is the one we use internally
- `default_value` -- the default value for the option
- `modify_dict` -- (optional, default true) indicates if `dict` should be modified
"""
function check_aliases_and_add_default_value!(dict::Dict{Symbol,Any}, dict_copy::Dict{Symbol,Any}, aliases::Vector{Symbol}, default_value::Any, modify_dict::Bool=true)
    check_aliases!(dict, dict_copy, aliases)

    if !haskey(dict_copy, aliases[1])
        # no alias and no value
        # TODO<notification> print message to user for each default option, depending on verbosity level
        if (modify_dict)
            dict[aliases[1]] = default_value
        end
        dict_copy[aliases[1]] = default_value
    end
end

"""
    check_valid_option_keywords(dict)

Validate that all input keywords are recognized.

### Input

- `dict` -- input dictionary
"""
function check_valid_option_keywords(dict)
    for kv_pair in dict
        key::Symbol = kv_pair.first
        if !in(key, available_keywords)
            error(get_unrecognized_key_message(key))
        end
    end
end

"""
    get_unrecognized_key_message(key)

Create an error message for an unrecognized option key.

### Input

- `key` -- unrecognized option key

### Output

The error message.
"""
function get_unrecognized_key_message(key)
    return "Unrecognized option '$key' found. See " *
        "`keys(Reachability.available_keywords.dict)` for all valid keywords."
end
