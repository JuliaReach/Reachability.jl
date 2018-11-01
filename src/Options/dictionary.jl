    import Base: merge, merge!, getindex, keys, haskey, values, setindex!, copy

    export Options, merge, merge!, getindex, haskey

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
