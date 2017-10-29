import Base: merge, getindex

export Options, merge, getindex

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

"""
Handling of the special trio `:δ`, `:N`, `:T`.
Usually two of them should be defined and the third one is automatically inferred.
If three of them are defined, they must be consistent.
If only `:T` is defined, we use `:N = 100`.

INPUT:

- ``dict``      -- a dictionary of options
- ``dict_copy`` -- a copy of the dictionary of options for internal names

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
    check_aliases!(dict, dict_copy, T_aliases)

    defined = 0
    if haskey(dict_copy, :δ)
        value = dict_copy[:δ]
        if !(value isa Float64 && value > 0.)
            error(string(value) * " is not a valid value for option $δ_aliases.")
        end
        defined += 1
    end
    if haskey(dict_copy, :N)
        value = dict_copy[:N]
        if !(value isa Int64 && value > 0)
            error(string(value) * " is not a valid value for option :N.")
        end
        defined += 1
    end
    if haskey(dict_copy, :T)
        value = dict_copy[:T]
        if !(value isa Float64 && value > 0.)
            error(string(value) * " is not a valid value for option $T_aliases.")
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
            error("Values for :δ, :N, and :T were defined, but they are inconsistent.")
        end
    else
        error("Need two of the following options: :δ, :N, :T")
    end
end


"""
    check_aliases!(dict, dict_copy, aliases)

This function has several purposes:

- translate aliases to the option that is used internally
- check that not several aliases were used at the same time

INPUT:

- ``dict`` -- a dictionary of options
- ``dict_copy`` -- a copy of the dictionary of options for internal names
- ``aliases`` -- option aliases; the first name is the one we use internally
- ``default_value`` -- the default value for the option
"""
function check_aliases!(dict::Dict{Symbol,Any}, dict_copy::Dict{Symbol,Any}, aliases::Vector{Symbol})
    # find aliases and check consistency in case several aliases have been used
    print_warning::Bool = false
    for alias in aliases
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
    check_aliases_and_add_default_value!(dict, dict_copy, aliases, default_value)

This function has several purposes:
- translate aliases to the option that is used internally
- check that not several aliases were used at the same time
- assign default values for undefined options

INPUT:

- ``dict`` -- a dictionary of options
- ``dict_copy`` -- a copy of the dictionary of options for internal names
- ``aliases`` -- option aliases; the first name is the one we use internally
- ``default_value`` -- the default value for the option
"""
function check_aliases_and_add_default_value!(dict::Dict{Symbol,Any}, dict_copy::Dict{Symbol,Any}, aliases::Vector{Symbol}, default_value::Any)
    check_aliases!(dict, dict_copy, aliases)

    if !haskey(dict, aliases[1])
        # no alias and no value
        # TODO<notification> print message to user for each default option, depending on verbosity level
        dict[aliases[1]] = default_value
        dict_copy[aliases[1]] = default_value
    end
end
