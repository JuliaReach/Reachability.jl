#============================================================
Struct with the dictionary of options and basic functionality
============================================================#

import Base: merge, merge!, getindex, keys, haskey, values, setindex!, copy,
             iterate, show

export Options,
       TwoLayerOptions, specified_keys, specified_values

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
    Options

Type that wraps two `Options` structs, one for specified options and one for
fallback defaults.

### Fields

- `specified` -- specified options
- `defaults`  -- default options

### Notes

It is possible to define `specified` options that are not contained in the
`defaults` options.

### Examples

```jldoctest
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
struct TwoLayerOptions
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

function getindex(𝑂::TwoLayerOptions, sym::Symbol)
    if haskey(𝑂.specified)
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
        print(io, "\n $key => $val")
    end
    print(io, "\nunspecified (default) options:")
    for (key, val) in 𝑂.defaults
        if !haskey(𝑂.specified, key)
            print(io, "\n $key => $val")
        end
    end
end
