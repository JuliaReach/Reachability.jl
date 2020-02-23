"""
    validate_solver_options_and_add_default_values!(options)

Validates that the given solver options are supported and adds default values for most
unspecified options.

### Input

- `options` -- an `Options` object, a dictionary of options

Supported options:

- `:verbosity`     -- controls logging output
- `:logfile`       -- name of a log file
- `:mode`          -- main analysis mode
- `:property`      -- `nothing`, a safety property, or a mapping from a location
                      to a safety property
- `:T`             -- time horizon; alias `:time_horizon`
- `:ε_proj`        -- error bound for the approximation of the states during
                      projection
- `:set_type_proj` -- set type for the approximation of the states during
                      projection
- `:coordinate_transformation` -- coordinate transformation method
- `:projection_matrix`         -- projection matrix
- `:project_reachset`          -- switch for applying projection
- `:max_jumps`     -- maximum number of discrete jumps in a hybrid automaton;
                      `-1` for deactivation
- `:fixpoint_check` -- check for a fixpoint when analyzing a hybrid automaton
- `:clustering`    -- clustering strategy when analyzing a hybrid automaton
- `:plot_vars`     -- variables for projection and plotting;
                      alias: `:output_variables`
- `:n`             -- system's dimension

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
    LOGGER = haskey(dict, :verbosity) ?
        configure_logger(dict[:verbosity]) :
        configure_logger()
    check_aliases_and_add_default_value!(dict, dict_copy, [:verbosity], getlevel(LOGGER))
    if haskey(dict, :logfile)
        add_file_logger(dict[:logfile])
        check_aliases!(dict, dict_copy, [:logfile])
    end

    # check for aliases and use default values for unspecified options
    check_aliases_and_add_default_value!(dict, dict_copy, [:mode], "reach")
    check_aliases_and_add_default_value!(dict, dict_copy, [:property], nothing)
    check_aliases_and_add_default_value!(dict, dict_copy, [:coordinate_transformation], "")
    check_aliases_and_add_default_value!(dict, dict_copy, [:projection_matrix], nothing)
    check_aliases_and_add_default_value!(dict, dict_copy, [:project_reachset], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:max_jumps], -1)
    check_aliases_and_add_default_value!(dict, dict_copy, [:clustering], :chull)
    check_aliases_and_add_default_value!(dict, dict_copy, [:fixpoint_check], :eager)
    check_aliases_and_add_default_value!(dict, dict_copy, [:n], nothing)
    check_aliases_and_add_default_value!(dict, dict_copy, [:transformation_matrix], nothing)
    check_aliases_and_add_default_value!(dict, dict_copy, [:plot_vars, :output_variables], [0, 1])
    check_aliases_and_add_default_value!(dict, dict_copy, [:ε_proj], Inf)
    check_aliases_and_add_default_value!(dict, dict_copy, [:set_type_proj], Hyperrectangle)

    # special option: T
    check_aliases!(dict, dict_copy, [:T, :time_horizon])

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
        elseif key == :logfile
            expected_type = String
        elseif key == :mode
            expected_type = String
            domain_constraints = (v::String  ->  v in ["reach", "check"])
            if value == "check" && dict_copy[:property] == nothing
                error("No property has been defined.")
            end
        elseif key == :property
            expected_type = Union{Nothing, Predicate, Dict{Int, <:Predicate}, Function, Dict{Int, <:Function}}
        elseif key == :T
            expected_type = Float64
            domain_constraints = (v::Float64  ->  v > 0.)
        elseif key == :ε_proj
            expected_type = Float64
            domain_constraints = (v  ->  v > 0.)
        elseif key == :set_type_proj
            expected_type = Union{Type{HPolygon}, Type{Hyperrectangle},
                                  Type{LazySets.Interval}}
        elseif key == :coordinate_transformation
            expected_type = String
            domain_constraints = (v::String  ->  v in ["", "schur"])
        elseif key == :projection_matrix
            expected_type = Union{AbstractMatrix, Nothing}
        elseif key == :project_reachset
            expected_type = Bool
        elseif key == :max_jumps
            expected_type = Int
            domain_constraints = (v::Int  ->  v >= -1)
        elseif key == :fixpoint_check
            expected_type = Symbol
            domain_constraints = (v::Symbol  ->  v in [:none, :eager, :lazy])
        elseif key == :clustering
            expected_type = Symbol
            domain_constraints = (v::Symbol  ->  v in [:chull, :none, :none_oa])
        elseif key == :plot_vars
            expected_type = Vector{Int}
            domain_constraints = (v::Vector{Int}  ->  length(v) == 2)
        elseif key == :n
            expected_type = Int
            domain_constraints = (v::Int  ->  v > 0)
        elseif key == :transformation_matrix
            expected_type = Any
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
    # ε-close approximation
    if dict_copy[:ε_proj] < Inf && dict_copy[:set_type_proj] != HPolygon
        throw(DomainError("ε-close approximation is only supported with the " *
                          "set type 'HPolygon'"))
    end

    return options_copy
end
