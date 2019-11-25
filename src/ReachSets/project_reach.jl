"""
    project_reach(Rsets, vars, options)

Projection of a reachability analysis result in 2D.

### Input

- `Rsets`   -- reachable states representation
- `vars`    -- variables to plot; two-dimensional index vector
- `options` -- options

### Notes

The `vars` argument is required even if the optional argument
`projection_matrix` is passed, because we also determine whether time is used as
a dimension from this variable.
"""
function project_reach(
        Rsets::Vector{<:AbstractReachSet{<:LazySets.LazySet{N}}},
        vars::Vector{Int64},
        options::AbstractOptions)::Vector{<:AbstractReachSet} where {N<:Real}
    # parse input
    @assert length(vars) == 2 "we only support projection to two dimensions"
    xaxis = vars[1]
    yaxis = vars[2]
    got_time = (xaxis == 0)
    if xaxis < 0
        throw(DomainError("value $xaxis for x variable not allowed"))
    elseif yaxis <= 0
        throw(DomainError("value $yaxis for y variable not allowed"))
    end
    n = options[:n]
    projection_matrix_high_dimensional = options[:projection_matrix]
    transformation_matrix = options[:transformation_matrix]

    # allocate output and define overapproximation function
    ε = options[:ε_proj]
    if ε < Inf
        oa = x -> overapproximate(x, HPolygon, ε)
        RsetsProj = Vector{ReachSet{HPolygon{N}}}(undef, length(Rsets))
    else
        set_type = options[:set_type_proj]
        oa = x -> overapproximate(x, set_type)
        RsetsProj = Vector{ReachSet{<:set_type{N}}}(undef, length(Rsets))
    end

    @inbounds for (i, rs) in enumerate(Rsets)
        X = set(rs)
        if projection_matrix_high_dimensional == nothing
            # use a simple projection to state variables
            if got_time
                # projection to a single state variable
                projection_matrix = sparse([1], [yaxis], [1.0], 1, n)
            else
                # projection to two state variables
                projection_matrix =
                    sparse([1, 2], [xaxis, yaxis], [1.0, 1.0], 2, n)
            end
        else
            @assert size(projection_matrix_high_dimensional, 1) == 1 "currently " *
                "we only support one-dimensional projection matrices"
            if got_time
                # create a 1-row matrix
                projection_matrix =
                    sparse(fill(1, n), 1:n, projection_matrix_high_dimensional[1, :], 1, n)
            else
                # create a 2-row matrix
                projection_matrix =
                    sparse(fill(2, n), 1:n, projection_matrix_high_dimensional[1, :], 2, n)
                projection_matrix[1, xaxis] = 1.0
            end
        end

        # apply optional transformation to projection matrix
        if transformation_matrix != nothing
            if got_time
                # add another dimension for time: block matrix [S 0; 0 1]
                transformation_matrix =
                    sparse(cat([1, 2], transformation_matrix, [1]))
            end
            projection_matrix = projection_matrix * transformation_matrix
        end

        # project set
        rs = project(rs, projection_matrix)
        projected = set(rs)

        # add time dimension
        t0 = time_start(rs)
        t1 = time_end(rs)
        if got_time
            time_interval = Interval(t0, t1)
            projected = CartesianProduct(time_interval, projected)
        end
        RsetsProj[i] = ReachSet(oa(projected), t0, t1)
    end

    return RsetsProj
end

"""
    project(Rsets, options)

Projects a sequence of sets according to the settings defined in the options.

### Input

- `Rsets`   -- solution of a reachability problem
- `options` -- options structure

### Notes

A projection matrix can be given in the options structure, or passed as a
dictionary entry.
"""
function project(Rsets::Vector{<:AbstractReachSet}, options::AbstractOptions)
    return project_reach(Rsets, options[:plot_vars], options)
end

project(reach_sol::ReachSolution) = project(reach_sol.Xk, reach_sol.options)

project(Rsets::Vector{<:AbstractReachSet}, options::Pair{Symbol,<:Any}...) =
    project(Rsets, Options(Dict{Symbol,Any}(options)))
