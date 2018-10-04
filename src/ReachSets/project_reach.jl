"""
    project_reach(plot_vars, n, δ, Rsets; [ε], [projection_matrix], [transformation_matrix])

Projection of a reachability analysis result in 2D.

### Input

- `plot_vars`         -- variables to plot; two-dimensional index vector
- `n`                 -- system dimension
- `δ`                 -- time discretization
- `Rsets`             -- reachable states representation
- `ε`                 -- (optional, default: `Inf`) error bound for the
                         approximation
- `set_type`          -- (optional, default: `Hyperrectangle`) set type for the
                         approximation
- `projection_matrix` -- (optional, default: `nothing`) projection matrix; if
                         not passed, the function computes `projection_matrix`
                         from `plot_vars`
- `transformation_matrix` -- (optional, default: `nothing`) transformation
                             matrix
- `output_function`   -- (optional, default: `false`) switch denoting whether
                         the passed set is one-dimensional, representing an
                         output function

### Notes

The `plot_vars` argument is required even if the optional argument
`projection_matrix` is passed, because we also determine whether time is used as
a dimension from this variable.
"""
function project_reach(plot_vars::Vector{Int64}, n::Int64, δ::Float64,
    Rsets::Vector{<:ReachSet{<:LazySets.CartesianProductArray{numeric_type}}};
    ε::Float64=Inf, set_type::Type{<:LazySet}=Hyperrectangle,
    projection_matrix::Union{AbstractMatrix, Void}=nothing,
    transformation_matrix::Union{AbstractMatrix, Void}=nothing,
    output_function::Bool=false
    )::Vector{<:ReachSet} where {numeric_type<:Real}

    # parse input
    assert(length(plot_vars) == 2)
    # first projection dimension
    xaxis = plot_vars[1]
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
    m = got_time ? n+1 : n
    if projection_matrix == nothing
        # projection to a state variable
        yaxis = plot_vars[2]
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
    if (transformation_matrix != nothing)
        if got_time
            # add another dimension for time: block matrix [S 0; 0 1]
            transformation_matrix = sparse(cat([1, 2], transformation_matrix, [1]))
        end
        projection_matrix = projection_matrix * transformation_matrix
    end

    N = length(Rsets)

    # allocate output and define overapproximation function
    if ε < Inf
        oa = x -> overapproximate(x, HPolygon, ε)
        RsetsProj = Vector{ReachSet{HPolygon{numeric_type}, numeric_type}}(N)
    else
        oa = x -> overapproximate(x, set_type)
        RsetsProj = Vector{ReachSet{set_type{numeric_type}, numeric_type}}(N)
    end

    if got_time
        radius = δ/2.0
        t = radius
        @inbounds for i in 1:N
            RsetsProj[i] = ReachSet(
                oa(projection_matrix *
                    CartesianProduct(Rsets[i].X, BallInf([t], radius))),
                Rsets[i].t_start, Rsets[i].t_end)
            t = t + δ
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
    project_reach(plot_vars, n, δ, Rsets; [ε], [projection_matrix],
                  [transformation_matrix])

This function projects a sequence of sets into the time variable, or can be
used to take a linear combination of the given variables.

The input `Rsets` is an array of sets (instead of a `CartesianProductArray`).
This array contains the collection of reach sets in 2D.

WARNING: 

It is assumed that the variable given in plot_vars belongs to the block computed
in the sequence of 2D sets `Rsets`.
"""
function project_reach(plot_vars::Vector{Int64}, n::Int64, δ::Float64,
    Rsets::Vector{<:ReachSet{<:LazySets.LazySet{numeric_type}}};
    ε::Float64=Inf, set_type::Type{<:LazySet}=Hyperrectangle,
    projection_matrix::Union{AbstractMatrix, Void}=nothing,
    transformation_matrix::Union{AbstractMatrix, Void}=nothing,
    output_function::Bool=false
    )::Vector{<:ReachSet} where {numeric_type<:Real}

    # parse input
    assert(length(plot_vars) == 2)
    # first projection dimension
    xaxis = plot_vars[1]
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
    if projection_matrix == nothing
        # projection to a state variable
        yaxis = plot_vars[2]
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
    if ε < Inf
        oa = x -> overapproximate(x, HPolygon, ε)
        RsetsProj = Vector{ReachSet{HPolygon{numeric_type}, numeric_type}}(N)
    else
        oa = x -> overapproximate(x, set_type)
        RsetsProj = Vector{ReachSet{set_type{numeric_type}, numeric_type}}(N)
    end

    if output_function
        radius = δ/2.0
        t = radius
        @inbounds for i in 1:N
            RsetsProj[i] = ReachSet(
                oa(CartesianProduct(BallInf([t], radius), Rsets[i].X)),
                Rsets[i].t_start, Rsets[i].t_end)
            t = t + δ
        end
    elseif got_time # x variable is 'time'
        radius = δ/2.0
        t = radius
        @inbounds for i in 1:N
            RsetsProj[i] = ReachSet(
                oa(projection_matrix *
                    CartesianProduct(Rsets[i].X, BallInf([t], radius))),
                Rsets[i].t_start, Rsets[i].t_end)
            t = t + δ
        end
    else
        @inbounds for i in 1:N
            RsetsProj[i] = oa(projection_matrix * Rsets[i].X)
        end
    end

    return RsetsProj
end
