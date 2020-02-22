#============================================================
Struct with the dictionary of options and basic functionality
============================================================#

import Base: merge, merge!, getindex, keys, haskey, values, setindex!, copy,
             iterate, show

export Options

"""
    AbstractOptions

Abstract supertype of options wrappers.
"""
abstract type AbstractOptions end

"""
    Options <: AbstractOptions

Type that wraps a dictionary used for options.

### Fields

- `dict` -- the wrapped dictionary
"""
struct Options <: AbstractOptions
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

```jldoctest
julia> op = Options(:T => 1.0)
Options(Dict{Symbol,Any}(:T => 1.0))

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

```jldoctest
julia> op = Options(:T => 1.0)
Options(Dict{Symbol,Any}(:T => 1.0))

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

```jldoctest
julia> op = Options(:T => 1.0)
Options(Dict{Symbol,Any}(:T => 1.0))

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
    iterate(op::Options)

Iterate over options.

### Input

- `op` -- options object
"""
iterate(op::Options) = iterate(op.dict)

iterate(op::Options, i::Int) = iterate(op.dict, i)

#=================================================
Struct for two-layered options with default values
=================================================#

"""
    TwoLayerOptions <: AbstractOptions

Type that wraps two `Options` structs, one for specified options and one for
fallback defaults.

### Fields

- `specified` -- specified options
- `defaults`  -- default options

### Notes

It is possible to define `specified` options that are not contained in the
`defaults` options.

### Examples

```julia
julia> def = Options(:o1 => "v1", :o2 => "v2");

julia> spec = Options(:o2 => "v2", :o3 => "v3");

julia> o = TwoLayerOptions(spec, def)
specified options:
 o2 => v2
 o3 => v3
unspecified (default) options:
 o1 => v1

```
"""
struct TwoLayerOptions <: AbstractOptions
    specified::Options
    defaults::Options
end

keys(𝑂::TwoLayerOptions) = keys(𝑂.defaults)

specified_keys(𝑂::TwoLayerOptions) = keys(𝑂.specified)

function values(𝑂::TwoLayerOptions)
    vals = values(𝑂.specified)
    for (key, val) in 𝑂.defaults
        if !haskey(𝑂.specified, key)
            push!(vals, val)
        end
    end
    return vals
end

specified_values(𝑂::TwoLayerOptions) = values(𝑂.specified)

function haskey(𝑂::TwoLayerOptions, key)
    return haskey(𝑂.specified, key) || haskey(𝑂.defaults, key)
end

function haskey_specified(𝑂::TwoLayerOptions, key)
    return haskey(𝑂.specified, key)
end

function getindex(𝑂::TwoLayerOptions, sym::Symbol)
    if haskey(𝑂.specified, sym)
        return getindex(𝑂.specified, sym)
    end
    return getindex(𝑂.defaults, sym)
end

function setindex!(𝑂::TwoLayerOptions, value, key)
    error("setting values in TwoLayerOptions is not allowed")
end

function show(io::IO, 𝑂::TwoLayerOptions)
    print(io, "specified options:")
    for (key, val) in 𝑂.specified
        if val == nothing
            val = "nothing"
        end
        print(io, "\n $key => $val")
    end
    print(io, "\nunspecified (default) options:")
    for (key, val) in 𝑂.defaults
        if !haskey(𝑂.specified, key)
            if val == nothing
                val = "nothing"
            end
            print(io, "\n $key => $val")
        end
    end
end

#====================================
Struct for specifying a single option
====================================#

"""
    OptionSpec{T}

Type that wraps the specification of an option.

### Fields

- `name`         -- name of the option as a symbol
- `default`      -- default value
- `aliases`      -- list of aliases (symbols)
- `domain_check` -- function for domain checks
- `info`         -- additional info text

### Examples

```julia
julia> os1 = OptionSpec(:option1, nothing)
option :option1 of type Any has default value 'nothing'

julia> os2 = OptionSpec(:option2, 0; aliases=[:othername, :yetanothername],
       domain=Number, domain_check=x->x>=0, info="the value is nonnegative")
option :option2 of type Number with aliases :othername and :yetanothername has default value '0' such that the value is nonnegative

```
"""
struct OptionSpec{T}
    name::Symbol
    default::T
    aliases::Vector{Symbol}
    domain_check::Function
    info::String

    OptionSpec(name::Symbol,
               default;
               aliases::Vector{Symbol}=Symbol[],
               domain::Type{T}=Any,
               domain_check::Function=(x -> true),
               info::String="") where {T} =
        new{T}(name, default, aliases, domain_check, info)
end

function show(io::IO, os::OptionSpec{T}) where {T}
    print(io, "option :", os.name, " of type ", string(T), " ")
    if !isempty(os.aliases)
        if length(os.aliases) == 1
            print(io, "with alias :", os.aliases[1], " ")
        elseif length(os.aliases) == 2
            print(io, "with aliases :", os.aliases[1], " and :", os.aliases[2],
                      " ")
        else
            print(io, "with aliases :")
            i = 1
            while i < length(os.aliases)
                print(io, os.aliases[i], ", :")
                i += 1
            end
            print(io, "and :", os.aliases[end], " ")
        end
    end
    print(io, "has default value '")
    if os.default == nothing
        print(io, "nothing")
    else
        print(io, os.default)
    end
    print(io, "'")
    if !isempty(os.info)
        print(io, "; ", os.info)
    end
end

function domain(::OptionSpec{T})::Type{T} where {T}
    return T
end

"""
    optionsspeclist_2_optionsspecmap(specs::AbstractVector{<:OptionSpec}
                                    )::Dict{Symbol, <:OptionSpec}

Crate a map from option names (including aliases) to the corresponding
specification.

### Input

- `specs` -- list of option specifications

### Output

A map from option names to option specification.

### Examples

```julia
julia> specs_list = [OptionSpec(:option1, nothing),
                     OptionSpec(:option2, nothing, aliases=[:op2, :op2_v2])];

julia> Reachability.optionsspeclist_2_optionsspecmap(specs_list)
Dict{Symbol,OptionSpec} with 4 entries:
  :option2 => option :option2 of type Any with aliases :op2 and :op2_v2 has default value 'nothing'
  :op2     => option :option2 of type Any with aliases :op2 and :op2_v2 has default value 'nothing'
  :op2_v2  => option :option2 of type Any with aliases :op2 and :op2_v2 has default value 'nothing'
  :option1 => option :option1 of type Any has default value 'nothing'
```
"""
function optionsspeclist_2_optionsspecmap(specs::AbstractVector{<:OptionSpec}
                                         )::Dict{Symbol, <:OptionSpec}
    function _add!(map, key, val)
        if haskey(map, key)
            throw(DomainError(":" * string(key), "duplicate option"))
        end
        map[key] = val
    end

    specs_map = Dict{Symbol, OptionSpec}()
    for spec in specs
        _add!(specs_map, spec.name, spec)
        for alias in spec.aliases
            _add!(specs_map, alias, spec)
        end
    end
    return specs_map
end

"""
    unify_aliases(𝑂::Options, specs::Dict{Symbol, <:OptionSpec})::Options

Crate a new options wrapper that uses the internal name for each option.

### Input

- `𝑂`     -- options wrapper
- `specs` -- list of option specifications

### Output

A new options wrapper.

### Examples

```jldoctest
julia> using Reachability: OptionSpec

julia> 𝑂 = Options(:option1 => 1.0, :op2 => "value");

julia> specs_list = [OptionSpec(:option1, nothing),
                     OptionSpec(:option2, nothing, aliases=[:op2, :op2_v2])];

julia> specs_map = Reachability.optionsspeclist_2_optionsspecmap(specs_list);

julia> Reachability.unify_aliases(𝑂, specs_map)
Options(Dict{Symbol,Any}(:option2 => "value",:option1 => 1.0))

```
"""
function unify_aliases(𝑂::Options, specs::Dict{Symbol, <:OptionSpec})::Options
    𝑂_new = Options()
    for (name, value) in 𝑂
        spec = get(specs, name, nothing)
        if spec == nothing
            warn("illegal option ':$name' detected; option will be ignored")
        else
            𝑂_new[spec.name] = value
        end
    end
    return 𝑂_new
end

"""
    validate_options(𝑂::Options, specs::Dict{Symbol, <:OptionSpec})

Validate the given options with respect to the given specifications.

### Input

- `𝑂`     -- options
- `specs` -- map from option names to option specifications

### Output

Nothing if the validation succeeds, and an error otherwise.
"""
function validate_options(𝑂::Options, specs::Dict{Symbol, <:OptionSpec})
    for (name, val) in 𝑂
        spec = specs[name]
        type = domain(spec)
        if !(val isa type)
            throw(DomainError(val, "option ':$name' must be of type $type"))
        end
        if !(spec.domain_check(val))
            throw(DomainError(val, "option ':$name' does not satisfy domain " *
                "constraint"))
        end
    end
    nothing
end

"""
    validate_and_wrap_options(𝑂::Options,
                              specs_list::AbstractVector{<:OptionSpec}
                             )::TwoLayerOptions

Normalize, validate, and merge given options with respect to default options.

### Input

- `𝑂`             -- options
- `specs`         -- list of option specifications
- `validation`    -- (optional, default: no-op) function for algorithm-specific
                     validation of the options
- `normalization` -- (optional, default: no-op) function for algorithm-specific
                     normalization of the options

### Output

A normalized options wrapper.

### Examples

```julia
julia> 𝑂 = Options(:option1 => 1.0, :op2 => "value");

julia> specs_list = [OptionSpec(:option1, nothing),
                     OptionSpec(:option2, nothing, aliases=[:op2, :op2_v2]),
                     OptionSpec(:option3, 2.0),];

julia> Reachability.validate_and_wrap_options(𝑂, specs_list)
specified options:
 option2 => value
 option1 => 1.0
unspecified (default) options:
 option3 => 2.0

```
"""
function validate_and_wrap_options(𝑂::Options,
                                   specs_list::AbstractVector{<:OptionSpec};
                                   validation::Function=(x -> nothing),
                                   normalization::Function=(x -> nothing)
                                  )::TwoLayerOptions
    specs_map = optionsspeclist_2_optionsspecmap(specs_list)
    𝑂_normalized = unify_aliases(𝑂, specs_map)
    validate_options(𝑂_normalized, specs_map)
    𝑂_default = Options(Dict{Symbol, Any}(
        spec.name => spec.default for spec in specs_list))
    𝑂_result = TwoLayerOptions(𝑂_normalized, 𝑂_default)
    normalization(𝑂_result)
    validation(𝑂_result)
    return 𝑂_result
end

"""
    print_option_spec(specs::AbstractVector{<:OptionSpec})

Pretty printing of a list of option specifications.

### Input

- `specs` -- list of option specifications

### Output

A string representation of the option specifications.
"""
function print_option_spec(specs::AbstractVector{<:OptionSpec})
    io = IOBuffer();
    for spec in specs
        println(io, spec)
    end
    return String(take!(io))
end
