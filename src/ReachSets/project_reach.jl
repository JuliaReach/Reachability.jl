"""
    project_reach(Rsets, vars, n, options)

Projection of a reachability analysis result in 2D.

### Input

- `Rsets`   -- reachable states representation
- `vars`    -- variables to plot; two-dimensional index vector
- `n`       -- system dimension
- `options` -- options

### Notes

The `vars` argument is required even if the optional argument
`projection_matrix` is passed, because we also determine whether time is used as
a dimension from this variable.
"""
function project_reach(
        Rsets::Vector{<:ReachSet{<:LazySets.LazySet{numeric_type}}},
        vars::Vector{Int64},
        n::Int64,
        options::AbstractOptions)::Vector{<:ReachSet} where {numeric_type<:Real}

    # parse input
    @assert(length(vars) == 2)
    if n == 2
        return project_2d_reach(Rsets, vars, n, options)
    end

    # first projection dimension
    xaxis = vars[1]
    if xaxis == 0
        got_time = true
        xaxis = n+1 # we add a new dimension for time
    else
        got_time = false
        if (xaxis <= 0 || xaxis > n)
            throw(DomainError())
        end
    end

    # build projection matrix
    projection_matrix = options[:projection_matrix]
    output_function = !options[:project_reachset]
    m = got_time ? n+1 : n
    if projection_matrix == nothing
        # projection to a state variable
        yaxis = vars[2]
        if (yaxis <= 0 || yaxis > n)
            throw(DomainError())
        end
        projection_matrix = sparse([1, 2], [xaxis, yaxis], [1.0, 1.0], 2, m)
    else
        @assert(size(projection_matrix) == (1,n))
        # make vector a 2-row matrix
        projection_matrix = sparse(fill(2, n), 1:n, projection_matrix[1, :], 2, m)
        projection_matrix[1, xaxis] = 1.0
    end

    # apply optional transformation to projection matrix
    if options[:transformation_matrix] != nothing
        transformation_matrix = options[:transformation_matrix]
        if got_time
            # add another dimension for time: block matrix [S 0; 0 1]
            transformation_matrix = sparse(cat([1, 2], transformation_matrix, [1]))
        end
        projection_matrix = projection_matrix * transformation_matrix
    end

    N = length(Rsets)

    # allocate output and define overapproximation function
    ε = options[:ε_proj]
    if ε < Inf
        oa = x -> overapproximate(x, HPolygon, ε)
        RsetsProj = Vector{ReachSet{HPolygon{numeric_type}, numeric_type}}(undef, N)
    else
        set_type = options[:set_type_proj]
        oa = x -> overapproximate(x, set_type)
        RsetsProj = Vector{ReachSet{set_type{numeric_type}, numeric_type}}(undef, N)
    end

    if got_time
        @inbounds for i in 1:N
            t0 = Rsets[i].t_start
            t1 = Rsets[i].t_end
            radius = (t1 - t0)/2.0
            RsetsProj[i] = ReachSet(
                oa(projection_matrix *
                    CartesianProduct(Rsets[i].X, BallInf([t0 + radius], radius))),
                Rsets[i].t_start, Rsets[i].t_end)
        end
    else
        @inbounds for i in 1:N
            RsetsProj[i] = ReachSet(
                oa(projection_matrix * Rsets[i].X),
                Rsets[i].t_start, Rsets[i].t_end)
        end
    end

    return RsetsProj
end

"""
    project_reach(Rsets, vars, n, options)

This function projects a sequence of sets into the time variable, or can be
used to take a linear combination of the given variables.

### Input

- `Rsets`   -- reachable states representation
- `vars`    -- variables to plot; two-dimensional index vector
- `n`       -- system dimension
- `options` -- options

### Notes

The input `Rsets` is an array of sets (instead of a `CartesianProductArray`).
This array contains the collection of reach sets in 2D.

It is assumed that the variable given in vars belongs to the block computed
in the sequence of 2D sets `Rsets`.
"""
function project_2d_reach(
        Rsets::Vector{<:ReachSet{<:LazySets.LazySet{numeric_type}}},
        vars::Vector{Int64},
        n::Int64,
        options::AbstractOptions)::Vector{<:ReachSet} where {numeric_type<:Real}

    # first projection dimension
    xaxis = vars[1]
    if xaxis == 0
        got_time = true
        xaxis = 3  # time is associated to dimension 3
    else
        got_time = false
        if (xaxis <= 0 || xaxis > n)
            throw(DomainError())
        end
        # map back to 2d block
        xaxis = iseven(xaxis) ? 2 : 1
    end

    # build projection matrix
    projection_matrix = options[:projection_matrix]
    output_function = !options[:project_reachset]
    if projection_matrix == nothing
        # projection to a state variable
        yaxis = vars[2]
        if (yaxis <= 0 || yaxis > n)
            throw(DomainError())
        end
        # map back to 2d block
        yaxis = iseven(yaxis) ? 2 : 1
        m = got_time ? 3 : 2
        projection_matrix = sparse([1, 2], [xaxis, yaxis], [1.0, 1.0], 2, m)
    elseif !output_function
        error("projection matrix not allowed for this algorithm")
    end

    N = length(Rsets)

    # allocate output and define overapproximation function
    ε = options[:ε_proj]
    if ε < Inf
        oa = x -> overapproximate(x, HPolygon, ε)
        RsetsProj = Vector{ReachSet{HPolygon{numeric_type}, numeric_type}}(undef, N)
    else
        set_type = options[:set_type_proj]
        oa = x -> overapproximate(x, set_type)
        RsetsProj = Vector{ReachSet{set_type{numeric_type}, numeric_type}}(undef, N)
    end

    if output_function
        @inbounds for i in 1:N
            t0 = Rsets[i].t_start
            t1 = Rsets[i].t_end
            radius = (t1 - t0)/2.0
            RsetsProj[i] = ReachSet(
                oa(CartesianProduct(BallInf([t0 + radius], radius), Rsets[i].X)),
                Rsets[i].t_start, Rsets[i].t_end)
        end
    elseif got_time # x variable is 'time'
        @inbounds for i in 1:N
            t0 = Rsets[i].t_start
            t1 = Rsets[i].t_end
            radius = (t1 - t0)/2.0
            RsetsProj[i] = ReachSet(
                oa(projection_matrix *
                    CartesianProduct(Rsets[i].X, BallInf([t0 + radius], radius))),
                Rsets[i].t_start, Rsets[i].t_end)
        end
    else
        @inbounds for i in 1:N
            RsetsProj[i] = ReachSet(oa(projection_matrix * Rsets[i].X),
                Rsets[i].t_start, Rsets[i].t_end, Rsets[i].k)
        end
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
function project(Rsets::Vector{<:ReachSet}, options::AbstractOptions)
    plot_vars = copy(options[:plot_vars])
    for i in 1:length(plot_vars)
        if plot_vars[i] != 0
            plot_vars[i] = options[:inout_map][plot_vars[i]]
        end
    end
    reduced_n = sum(x -> x != 0, options[:inout_map])
    output_function = !options[:project_reachset]
    RsetsProj = project_reach(Rsets, plot_vars, reduced_n, options)
end

project(reach_sol::ReachSolution) = project(reach_sol.Xk, reach_sol.options)

project(Rsets::Vector{<:ReachSet}, options::Pair{Symbol,<:Any}...) =
    project(Rsets, Options(Dict{Symbol,Any}(options)))
