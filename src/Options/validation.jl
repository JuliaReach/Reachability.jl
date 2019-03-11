#=============================================
Validation methods for dictionaries of options
==============================================#

available_keywords = Set{Symbol}([])

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
