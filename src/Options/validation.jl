#=============================================
Validation methods for dictionaries of options
==============================================#

available_keywords = Set{Symbol}([])

"""
    check_and_add_δ_N_T!(dict::Dict{Symbol,Any}, dict_copy::Dict{Symbol,Any})

Handling of the special trio `:δ`, `:N`, `:T`.
Usually two of them should be defined and the third one is automatically inferred.
If three of them are defined, they must be consistent.
If only `:T` is defined, we use `:N = 100`.

### Input

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
    check_and_add_plot_vars!(dict::Dict{Symbol,Any},
                             dict_copy::Dict{Symbol,Any})

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
