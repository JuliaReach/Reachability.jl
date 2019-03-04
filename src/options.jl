import Base: merge, merge!, getindex, keys, haskey, values, setindex!, copy

export Options, merge, merge!, getindex, haskey

available_keywords = Set{Symbol}([])

"""
    Options

Type that wraps a dictionary used for options.

### Fields

- `dict` -- the wrapped dictionary
"""
struct Options
    dict::Dict{Symbol,Any}
    Options(args::Pair{Symbol,<:Any}...) = new(Dict{Symbol,Any}(args))
    Options(dict::Dict{Symbol,Any}) = new(dict)
end

"""
    keys(op::Options)

Return the keys of the given options object.

### Input

- `op`    -- options object

### Examples

Obtain the keys of some options with one element:

```jldoctest options_setindex
julia> op = Options(:T=>1.0)
Reachability.Options(Dict{Symbol,Any}(Pair{Symbol,Any}(:T, 1.0)))

julia> collect(keys(op))
1-element Array{Symbol,1}:
 :T
```
"""
keys(op::Options) = keys(op.dict)

"""
    values(op::Options)

Return the values of the given options object.

### Input

- `op`    -- options object

### Examples

Obtain the values of some options with one element:

```jldoctest options_setindex
julia> op = Options(:T=>1.0)
Reachability.Options(Dict{Symbol,Any}(Pair{Symbol,Any}(:T, 1.0)))

julia> collect(values(op))
1-element Array{Any,1}:
 1.0
```
"""
values(op::Options) = values(op.dict)

"""
    merge(op1, opn)

Merges two `Options` objects by just falling back to the wrapped `Dict` fields.
Values are inserted in the order in which the function arguments occur, i.e.,
for conflicting keys a later object overrides a previous value.

### Input

- `op1` -- first options object
- `opn` -- list of options objects

### Output

An `Options` object.

### Notes

This function makes a copy of the dictionary and does not modify its first
argument.
"""
function merge(op1::Options, opn::Options...)::Options
    merge!(copy(op1), opn...)
end

"""
    merge!(op1, opn)

Updates the first argument options `op1` with options from `opn`.

### Input

- `op1` -- first options object
- `opn` -- list of options objects

### Output

An `Options` object.

### Algorithm

Merges two `Options` objects by just falling back to the wrapped `Dict` fields.
Values are inserted in the order in which the function arguments occur, i.e.,
for conflicting keys a later object overrides a previous value.
"""
function merge!(op1::Options, opn::Options...)::Options
    dict = op1.dict
    for i in 1 : length(opn)
        merge!(dict, opn[i].dict)
    end
    return Options(dict)
end

"""
    getindex(op::Options, sym::Symbol)

Returns the value stored for key `sym`.

### Input

- `op`  -- options object
- `sym` -- key
"""
function getindex(op::Options, sym::Symbol)
    return getindex(op.dict, sym)
end

"""
    setindex!(op, value, key)

Store the given value at the given key in the options.

### Input

- `op`    -- options object
- `value` -- value
- `key`   -- key

### Examples

Create an empty options object and add an input:

```jldoctest options_setindex
julia> Options()
Reachability.Options(Dict{Symbol,Any}())

julia> op[:T] = 1.0
1.0
```
"""
setindex!(op::Options, value, key) = setindex!(op.dict, value, key)

"""
    copy(op)

Create a shallow copy of the given options.

### Input

- `op`    -- options object

### Output

A new `Options` instance whose dictionary is a copy of `op`'s dictionary.
"""
copy(op::Options) = Options(copy(op.dict))

"""
    haskey(op::Options, key)

Determine whether the given options has a mapping for a given key.

### Input

- `op` -- options object

### Output

`true` if `op` contains the option `key` and `false` otherwise.
"""
haskey(op::Options, key) = haskey(op.dict, key)

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
- `:approx_model`  -- model for bloating/continuous time analysis
- `:property`      -- a safety property
- `:δ`             -- time step; alias: `:sampling_time`
- `:N`             -- number of time steps
- `:T`             -- time horizon; alias `:time_horizon`
- `:algorithm`     -- algorithm backend
- `:vars`          -- variables of interest
- `:global_vars`          -- variables of interest for HS
- `:partition`     -- block partition; elements are 2D vectors containing the
                      start and end of a block
- `:block_types`   -- short hand to set `:block_types_init` and
                      `:block_types_iter`
- `:block_types_init` -- set type for the approximation of the initial states
                         for each block
- `:block_types_iter` -- set type for the approximation of the states ``X_k``,
                         ``k>0``, for each block
- `:ε`             -- short hand to set `:ε_init` and `:ε_iter`
- `:set_type`      -- short hand to set `:set_type_init` and `:set_type_iter`
- `:ε_init`        -- error bound for the approximation of the initial states
                      (during decomposition)
- `:set_type_init` -- set type for the approximation of the initial states
                      (during decomposition)
- `:ε_iter`        -- error bound for the approximation of the states ``X_k``,
                      ``k>0``
- `:set_type_iter` -- set type for the approximation of the states ``X_k``,
                      ``k>0``
- `:ε_proj`        -- error bound for the approximation of the states during
                      projection
- `:set_type_proj` -- set type for the approximation of the states during
                      projection
- `:lazy_expm`     -- switch for using lazy matrix exponential all the time
- `:assume_sparse` -- switch for sparse matrices
- `:pade_expm`     -- switch for using Pade approximant method
- `:lazy_X0`       -- switch for keeping the initial states a lazy set
- `:lazy_sih`      -- switch for using a lazy symmetric interval hull during the
                      discretization
- `:template_directions`       -- short hand to set `template_directions_init`
                                  and `template_directions_iter`
- `:template_directions_init`  -- directions to use for the approximation of the
                                  initial states (during decomposition)
- `:template_directions_iter`  -- directions to use for the approximation of the
                                  states ``X_k``, ``k>0``, for each block
- `:coordinate_transformation` -- coordinate transformation method
- `:assume_homogeneous`        -- switch for ignoring inputs
- `:projection_matrix`         -- projection matrix
- `:project_reachset`          -- switch for applying projection
- `:eager_checking`            -- switch for early terminating property checks
- `:lazy_inputs_interval`      -- length of interval in which the inputs are
                                  handled as a lazy set (``-1`` for 'never');
                                  generally may also be a predicate over indices
- `:lazy_expm_discretize`      -- switch to use lazy matrix exponential in the
                                  discretization phase (see also `:lazy_expm`)
- `:max_jumps`     -- maximum number of discrete jumps in a hybrid automaton
- `:fixpoint_check` -- check for a fixpoint when analyzing a hybrid automaton
- `:clustering`    -- clustering strategy when analyzing a hybrid automaton
- `:plot_vars`     -- variables for projection and plotting;
                      alias: `:output_variables`
- `:n`             -- system's dimension
- `:steps`         -- steps for which we need to compute continuous post

Internal options (inputs are ignored or even illegal):

- `:blocks`        -- list of all interesting block indices in the partition

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
    check_aliases_and_add_default_value!(dict, dict_copy, [:approx_model], "forward")
    check_aliases_and_add_default_value!(dict, dict_copy, [:property], nothing)
    check_aliases_and_add_default_value!(dict, dict_copy, [:algorithm], "explicit")
    check_aliases_and_add_default_value!(dict, dict_copy, [:vars], 1:options.dict[:n])
    check_aliases_and_add_default_value!(dict, dict_copy, [:global_vars], Vector{Int}())
    check_aliases_and_add_default_value!(dict, dict_copy, [:lazy_expm], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:assume_sparse], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:pade_expm], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:lazy_X0], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:lazy_sih], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:coordinate_transformation], "")
    check_aliases_and_add_default_value!(dict, dict_copy, [:assume_homogeneous], false)
    check_aliases_and_add_default_value!(dict, dict_copy, [:projection_matrix], nothing)
    check_aliases_and_add_default_value!(dict, dict_copy, [:project_reachset], dict_copy[:projection_matrix] == nothing)
    check_aliases_and_add_default_value!(dict, dict_copy, [:eager_checking], true)
    check_aliases_and_add_default_value!(dict, dict_copy, [:lazy_expm_discretize],
                                         dict_copy[:lazy_expm])
    check_aliases_and_add_default_value!(dict, dict_copy, [:max_jumps], 5)
    check_aliases_and_add_default_value!(dict, dict_copy, [:clustering], :chull)
    check_aliases_and_add_default_value!(dict, dict_copy, [:fixpoint_check], :eager)
    check_aliases_and_add_default_value!(dict, dict_copy, [:n], nothing)
    check_aliases_and_add_default_value!(dict, dict_copy, [:steps], Vector{Int}())

    # special options: δ, N, T
    check_and_add_δ_N_T!(dict, dict_copy)

    # special option: plot_vars
    check_and_add_plot_vars!(dict, dict_copy)

    # special options: ε, ε_init, ε_iter, ε_proj,
    #                  set_type, set_type_init, set_type_iter, set_type_proj,
    #                  template_directions, template_directions_init,
    #                  template_directions_iter
    check_and_add_approximation!(dict, dict_copy)

    # special options: partition, block_types, block_types_init, block_types_iter
    check_and_add_partition_block_types!(dict, dict_copy)

    # special option: lazy_inputs_interval
    check_and_add_lazy_inputs_interval!(dict, dict_copy)

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
        elseif key == :approx_model
            expected_type = String
            domain_constraints = (v::String  ->  v in ["forward", "backward",
                                                       "firstorder", "nobloating"])
        elseif key == :property
            expected_type = Union{Property, Nothing}
        elseif key == :δ
            expected_type = Float64
            domain_constraints = (v::Float64  ->  v > 0.)
        elseif key == :N
            expected_type = Int
            domain_constraints = (v::Int  ->  v > 0)
        elseif key == :T
            expected_type = Float64
            domain_constraints = (v::Float64  ->  v > 0.)
        elseif key == :algorithm
            expected_type = String
            domain_constraints = (v::String  ->  v in ["explicit", "wrap"])
        elseif key == :vars
            expected_type = AbstractVector{Int}
            domain_constraints =
                (v::AbstractVector{Int}  ->  length(v) > 0 && all(x -> x > 0, v))
        elseif key == :global_vars
            expected_type = AbstractVector{Int}
            domain_constraints =
                (v::AbstractVector{Int}  ->  length(v) >= 0 && all(x -> x > 0, v))
        elseif key == :partition
            expected_type = AbstractVector{<:AbstractVector{Int}}
        elseif key == :block_types
            expected_type = Union{Nothing,
                Dict{Type{<:LazySet}, AbstractVector{<:AbstractVector{Int}}}}
        elseif key == :block_types_init
            expected_type =
                Dict{Type{<:LazySet}, AbstractVector{<:AbstractVector{Int}}}
        elseif key == :block_types_iter
            expected_type =
                Dict{Type{<:LazySet}, AbstractVector{<:AbstractVector{Int}}}
        elseif key == :blocks
            expected_type = AbstractVector{Int}
            domain_constraints =
                (v::AbstractVector{Int}  ->  length(v) > 0 && all(x -> x > 0, v))
        elseif key == :ε
            expected_type = Float64
            domain_constraints = (v  ->  v > 0.)
        elseif key == :ε_init
            expected_type = Float64
            domain_constraints = (v  ->  v > 0.)
        elseif key == :ε_iter
            expected_type = Float64
            domain_constraints = (v  ->  v > 0.)
        elseif key == :ε_proj
            expected_type = Float64
            domain_constraints = (v  ->  v > 0.)
        elseif key == :set_type
            expected_type = Union{Type{HPolygon}, Type{Hyperrectangle},
                                  Type{LazySets.Interval}}
        elseif key == :set_type_init
            expected_type = Union{Type{HPolygon}, Type{Hyperrectangle},
                                  Type{LazySets.Interval}}
        elseif key == :set_type_iter
            expected_type = Union{Type{HPolygon}, Type{Hyperrectangle},
                                  Type{LazySets.Interval}}
        elseif key == :set_type_proj
            expected_type = Union{Type{HPolygon}, Type{Hyperrectangle},
                                  Type{LazySets.Interval}}
        elseif key == :lazy_expm
            expected_type = Bool
        elseif key == :assume_sparse
            expected_type = Bool
        elseif key == :pade_expm
            expected_type = Bool
        elseif key == :lazy_X0
            expected_type = Bool
        elseif key == :lazy_sih
            expected_type = Bool
        elseif key == :template_directions
            expected_type = Symbol
            domain_constraints = (v::Symbol  ->  v in [:box, :oct, :boxdiag,
                                                       :nothing])
        elseif key == :template_directions_init
            expected_type = Symbol
            domain_constraints = (v::Symbol  ->  v in [:box, :oct, :boxdiag,
                                                       :nothing])
        elseif key == :template_directions_iter
            expected_type = Symbol
            domain_constraints = (v::Symbol  ->  v in [:box, :oct, :boxdiag,
                                                       :nothing])
        elseif key == :coordinate_transformation
            expected_type = String
            domain_constraints = (v::String  ->  v in ["", "schur"])
        elseif key == :assume_homogeneous
            expected_type = Bool
        elseif key == :projection_matrix
            expected_type = Union{AbstractMatrix, Nothing}
        elseif key == :project_reachset
            expected_type = Bool
        elseif key == :eager_checking
            expected_type = Bool
        elseif key == :lazy_inputs_interval
            expected_type = Union{Int, Function, Nothing}
            domain_constraints = (v  ->  !(v isa Int) || v >= -1)
        elseif key == :lazy_expm_discretize
            expected_type = Bool
            if !value && dict_copy[:lazy_expm]
                error("cannot use option $(:lazy_expm) with deactivated " *
                      "option $(:lazy_expm_discretize)")
            end
        elseif key == :max_jumps
            expected_type = Int
            domain_constraints = (v::Int  ->  v >= 0)
        elseif key == :fixpoint_check
            expected_type = Symbol
            domain_constraints = (v::Symbol  ->  v in [:none, :eager, :lazy])
        elseif key == :clustering
            expected_type = Symbol
            domain_constraints = (v::Symbol  ->  v in [:chull, :none])
        elseif key == :plot_vars
            expected_type = Vector{Int}
            domain_constraints = (v::Vector{Int}  ->  length(v) == 2)
        elseif key == :n
            expected_type = Int
            domain_constraints = (v::Int  ->  v > 0)
        elseif key == :steps
            expected_type = AbstractVector{Int}
            domain_constraints =
                (v::AbstractVector{Int}  ->  length(v) >= 0 && all(x -> x > 0, v))
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
        if !(value isa Int && value > 0)
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
            N = (Int)(ceil(dict_copy[:T] / dict_copy[:δ]))
            dict[:N] = N
            dict_copy[:N] = N
        end
        if !haskey(dict_copy, :T)
            T = dict_copy[:N] * dict_copy[:δ]
            dict[:T] = T
            dict_copy[:T] = T
        end
    elseif defined == 3
        N_computed = (Int)(ceil(dict_copy[:T] / dict_copy[:δ]))
        if N_computed != dict_copy[:N]
            error("Values for :δ ($(dict_copy[:δ])), :N ($(dict_copy[:N])), and :T ($(dict_copy[:T])) were defined, but they are inconsistent.")
        end
    else
        error("Need two of the following options: :δ, :N, :T (or at least :T and using default values)")
    end
end

"""
Handling of the special option `:plot_vars`.

### Input

- `dict`      -- dictionary of options
- `dict_copy` -- copy of the dictionary of options for internal names

### Notes:

If no value is given, we take the first two dimensions from `:vars`.
If `:vars` has only one element, we use time for the other.
"""
function check_and_add_plot_vars!(dict::Dict{Symbol,Any},
                                 dict_copy::Dict{Symbol,Any})
    check_aliases!(dict, dict_copy, [:plot_vars, :output_variables])
    if !haskey(dict_copy, :plot_vars)
        vars = dict_copy[:vars]
        if length(vars) == 1
            plot_vars = [0, vars[1]]
        else
            plot_vars = [vars[1], vars[2]]
        end
        dict_copy[:plot_vars] = plot_vars
    else
        # sanity check
        vars = dict_copy[:vars]
        for i in dict_copy[:plot_vars]
            if i != 0 && i ∉ vars
                error("$(dict_copy[:plot_vars]) is not a subset of the variables $vars")
            end
        end
    end
end

"""
    check_and_add_approximation!(dict::Dict{Symbol,Any},
                                 dict_copy::Dict{Symbol,Any})

Handling of the special options `:ε` and `:set_type` resp. the subcases
`:ε_init`, `:ε_iter`, `:ε_proj`, `:set_type_init`, `:set_type_iter`,
`:set_type_proj`, `:template_directions`, `:template_directions_init`, and
`:template_directions_iter`.

### Input

- `dict`      -- dictionary of options
- `dict_copy` -- copy of the dictionary of options for internal names

### Notes:

`:ε` and `:set_type`, if defined, are used as fallbacks (if the more specific
options are undefined).
If both `:ε_init` and `:set_type_init` are defined, they must be consistent.
"""
function check_and_add_approximation!(dict::Dict{Symbol,Any},
                                      dict_copy::Dict{Symbol,Any})
    # fallback options
    check_aliases!(dict, dict_copy, [:ε])
    check_aliases!(dict, dict_copy, [:set_type])
    check_aliases!(dict, dict_copy, [:template_directions])
    ε = haskey(dict, :ε) ? dict[:ε] : Inf

    # fallback set type
    if haskey(dict, :set_type)
        # use the provided set type
        set_type = dict[:set_type]
    elseif ε < Inf
        # use polygons
        set_type = HPolygon
    elseif !haskey(dict, :partition)
        # use intervals
        set_type = Interval
        set_type_proj = Interval
    else
        # use hyperrectangles
        set_type = Hyperrectangle
    end

    ε_init = (haskey(dict, :set_type_init) && dict[:set_type_init] == HPolygon) ||
             (!haskey(dict, :set_type_init) && set_type == HPolygon) ? ε : Inf
    check_aliases_and_add_default_value!(dict, dict_copy, [:ε_init], ε_init)

    set_type_init = dict_copy[:ε_init] < Inf ? HPolygon : set_type
    check_aliases_and_add_default_value!(dict, dict_copy, [:set_type_init], set_type_init)

    template_directions_init = haskey(dict, :template_directions_init) ?
        dict[:template_directions_init] : haskey(dict, :template_directions) ?
        dict[:template_directions] : :nothing
    check_aliases_and_add_default_value!(dict, dict_copy,
        [:template_directions_init], template_directions_init)

    ε_iter = (haskey(dict, :set_type_iter) && dict[:set_type_iter] == HPolygon) ||
             (!haskey(dict, :set_type_iter) && set_type == HPolygon) ? ε : Inf
    check_aliases_and_add_default_value!(dict, dict_copy, [:ε_iter], ε_iter)

    set_type_iter = dict_copy[:ε_iter] < Inf ? HPolygon : set_type
    check_aliases_and_add_default_value!(dict, dict_copy, [:set_type_iter], set_type_iter)

    template_directions_iter = haskey(dict, :template_directions_iter) ?
        dict[:template_directions_iter] : haskey(dict, :template_directions) ?
        dict[:template_directions] : :nothing
    check_aliases_and_add_default_value!(dict, dict_copy,
        [:template_directions_iter], template_directions_iter)

    ε_proj = (haskey(dict, :set_type_proj) && dict[:set_type_proj] == HPolygon) ||
             (!haskey(dict, :set_type_proj) && set_type == HPolygon) ? ε :
             min(dict_copy[:ε_init], dict_copy[:ε_iter], ε)
    check_aliases_and_add_default_value!(dict, dict_copy, [:ε_proj], ε_proj)

    set_type_proj = dict_copy[:ε_proj] < Inf ? HPolygon :
        set_type == Interval && length(dict_copy[:plot_vars]) > 1 ?
        Hyperrectangle : set_type
    check_aliases_and_add_default_value!(dict, dict_copy, [:set_type_proj], set_type_proj)

    @assert (dict_copy[:ε_init] == Inf || dict_copy[:set_type_init] == HPolygon) &&
        (dict_copy[:ε_iter] == Inf || dict_copy[:set_type_iter] == HPolygon) &&
        (dict_copy[:ε_proj] == Inf || dict_copy[:set_type_proj] == HPolygon) (
            "ε-close approximation is only supported with the HPolygon set type")
end

"""
    check_and_add_partition_block_types!(dict::Dict{Symbol,Any},
                                         dict_copy::Dict{Symbol,Any})

Handling of the special options `:partition` and `:block_types` resp. the
subcases `:block_types_init` and `:block_types_iter`.

### Input

- `dict`      -- dictionary of options
- `dict_copy` -- copy of the dictionary of options for internal names
"""
function check_and_add_partition_block_types!(dict::Dict{Symbol,Any},
                                              dict_copy::Dict{Symbol,Any})
    check_aliases!(dict, dict_copy, [:partition])
    if !haskey(dict_copy, :partition)
        partition = [[i] for i in 1:dict[:n]]
        dict_copy[:partition] =  partition
    end

    block_types = nothing
    if haskey(dict, :block_types)
        for (key, value) in dict[:block_types]
            @assert key <: LazySet "the keys of the `partition` dictionary should be lazy sets"
            @assert typeof(value) <: AbstractVector{<:AbstractVector{Int}} "the keys of the `partition` dictionary should be vectors of vectors"
        end
        block_types = convert(Dict{Type{<:LazySet}, AbstractVector{<:AbstractVector{Int}}}, dict[:block_types])
    elseif haskey(dict, :set_type)
        block_types = Dict{Type{<:LazySet}, AbstractVector{<:AbstractVector{Int}}}(dict[:set_type] => copy(dict_copy[:partition]))
    end
    check_aliases_and_add_default_value!(dict, dict_copy, [:block_types], block_types)
    dict_copy[:block_types] = block_types

    block_types_init = haskey(dict, :block_types_init) ?
        dict[:block_types_init] :
        block_types != nothing ? block_types :
            Dict{Type{<:LazySet}, AbstractVector{<:AbstractVector{Int}}}(
                dict_copy[:set_type_init] => copy(dict_copy[:partition])
            )
    check_aliases_and_add_default_value!(dict, dict_copy, [:block_types_init],
                                         block_types_init)

    block_types_iter = haskey(dict, :block_types_iter) ?
        dict[:block_types_iter] :
        block_types != nothing ? block_types :
            Dict{Type{<:LazySet}, AbstractVector{<:AbstractVector{Int}}}(
                dict_copy[:set_type_iter] => copy(dict_copy[:partition])
            )
    check_aliases_and_add_default_value!(dict, dict_copy, [:block_types_iter],
                                         block_types_iter)

    # TODO add check that arguments are consistent

    # compute :blocks
    if haskey(dict, :blocks)
        error("option :blocks is not allowed")
    end
    dict_copy[:blocks] = compute_blocks(dict_copy[:vars], dict_copy[:partition])
end

function compute_blocks(vars, partition)
    blocks = Vector{Int}()
    sizehint!(blocks, length(vars))
    next = 0
    var_idx = 1
    for (i, block) in enumerate(partition)
        next += length(block)
        if vars[var_idx] <= next
            push!(blocks, i)
            var_idx += 1
            while var_idx <= length(vars) && vars[var_idx] <= next
                var_idx += 1
            end
            if var_idx > length(vars)
                break
            end
        end
    end
    @assert var_idx == length(vars) + 1
    sizehint!(blocks, length(blocks))
    return blocks
end

"""
    check_and_add_lazy_inputs_interval!(dict::Dict{Symbol,Any},
                                        dict_copy::Dict{Symbol,Any})

Handling of the special option `:lazy_inputs_interval`.

### Input

- `dict`      -- dictionary of options
- `dict_copy` -- copy of the dictionary of options for internal names

### Notes

Originally, this option was used to overapproximate a lazy set of inputs at
every multiple of ``k`` steps, where ``k`` was the controllable option; hence
the term `interval` in the name.
Meanwhile, this option was generalized to a predicate over integers.
The predicate returns `true` iff the lazy set should be overapproximated.
However, we still support index numbers as input, which are then translated into
a predicate.

The default value for this option is `nothing` or the index ``0``, which results
in overapproximation in each iteration.
Note that internally this has a different effect than passing the predicate that
always returns `true` because functions cannot be checked for equality.
Thus this option should simply not be set when not wanting this behavior (or the
value ``0`` should be used).

The input ``-1`` is interpreted as the predicate that always returns `false`.
"""
function check_and_add_lazy_inputs_interval!(dict::Dict{Symbol,Any},
                                             dict_copy::Dict{Symbol,Any})
    check_aliases!(dict, dict_copy, [:lazy_inputs_interval])
    if haskey(dict_copy, :lazy_inputs_interval)
        if dict_copy[:lazy_inputs_interval] == -1
            dict_copy[:lazy_inputs_interval] = (k -> false)
        elseif dict_copy[:lazy_inputs_interval] == 0
            dict_copy[:lazy_inputs_interval] = nothing
        elseif dict_copy[:lazy_inputs_interval] isa Int
            m = dict_copy[:lazy_inputs_interval]
            dict_copy[:lazy_inputs_interval] = (k -> k % m == 0)
        end
    else
        dict_copy[:lazy_inputs_interval] = nothing
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

### Input

- `dict`          -- a dictionary of options
- `dict_copy`     -- a copy of the dictionary of options for internal names
- `aliases`       -- option aliases; the first name is the one we use internally
- `default_value` -- the default value for the option
- `modify_dict`   -- (optional, default: `false`) indicates if `dict` should be
                     modified
"""
function check_aliases_and_add_default_value!(dict::Dict{Symbol,Any}, dict_copy::Dict{Symbol,Any}, aliases::Vector{Symbol}, default_value::Any, modify_dict::Bool=false)
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
