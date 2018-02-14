"""
    project_reach(plot_vars, n, δ, Rsets, [algorithm]; [ɛ], [projection_matrix], [transformation_matrix])

Projection of a reachability analysis result in 2D.

INPUT:

- ``plot_vars``   -- variables to plot; two-dimensional index array
- ``n``           -- system dimension
- ``δ``           -- time discretization
- ``Rsets``       -- reachable states representation
- ``algorithm``   -- (optional, default: ``'explicit'``) reachability algorithm
                     backend, see ``available_algorithms``
- ``ɛ``           -- (optional, default: Inf) error tolerance (possibly ignored
                     for some `set_type` arguments)
- ``set_type``    -- (optional, default: `HPolygon`) type of set that is used
                     for overapproximation in 2D
- ``projection_matrix`` -- (optional, default: ``nothing``) projection matrix;
                     if not passed, the function computes `projection_matrix`
                     from plot_vars
- ``transformation_matrix`` -- (optional, default: ``nothing``) transformation
                               matrix

NOTES:

The ``plot_vars`` argument is required even if the optional argument
``projection_matrix`` is passed, because we also determine whether ``'time'`` is
used as a dimension from this variable.
"""
function project_reach(plot_vars::Vector{Int64}, n::Int64,
    δ::Float64, Rsets::Vector{<:LazySets.CartesianProductArray},
    algorithm::String="explicit";
    ɛ::Float64=Inf, set_type::Type=HPolygon,
    projection_matrix::Union{SparseMatrixCSC{Float64,Int64}, Void}=nothing,
    transformation_matrix::Union{SparseMatrixCSC{Float64,Int64}, Void}=nothing
    )::Vector{<:set_type}

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

    # allocate output
    RsetsProj = Vector{set_type}(N)
    oa_arg = set_type == HPolygon ? ɛ : set_type

    if got_time
        radius = δ/2.0
        t = radius
        @inbounds for i in 1:N
            RsetsProj[i] = overapproximate(projection_matrix *
                CartesianProduct(Rsets[i], BallInf([t], radius)), oa_arg)
            t = t + δ
        end
    else
        @inbounds for i in 1:N
            RsetsProj[i] = overapproximate(projection_matrix * Rsets[i], oa_arg)
        end
    end

    return RsetsProj
end

"""
    project_reach(plot_vars, n, δ, Rsets, [algorithm]; [ɛ], [projection_matrix], [transformation_matrix])

This algorithm projects a sequence of polygons into the time variable, or can be
used to take a linear combination of the given variables.

The input Rsets is an array of polygons (instead of a CartesianProductArray).
This array of polygons contains the collection of reachable sets in 2d. 

WARNING: 

It is assumed that the variable given in plot_vars belongs to the block computed
in the sequence of 2d polygons, Rsets.
"""
function project_reach(plot_vars::Vector{Int64}, n::Int64, δ::Float64,
    Rsets::Vector{<:LazySets.LazySet}, algorithm::String;
    ɛ::Float64=Inf, set_type::Type=HPolygon,
    projection_matrix::Union{SparseMatrixCSC{Float64,Int64}, Void}=nothing,
    transformation_matrix::Union{SparseMatrixCSC{Float64,Int64}, Void}=nothing
    )::Vector{<:set_type}

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
    else
        error("projection matrix not allowed for this algorithm")
    end

    N = length(Rsets)

    # allocate output
    RsetsProj = Vector{set_type}(N)
    oa_arg = set_type == HPolygon ? ɛ : set_type

    if got_time # x variable is 'time'
        radius = δ/2.0
        t = radius
        @inbounds for i in 1:N
            RsetsProj[i] = overapproximate(projection_matrix * 
                CartesianProduct(Rsets[i], BallInf([t], radius)), oa_arg)
            t = t + δ
        end
    else
        @inbounds for i in 1:N
            RsetsProj[i] = overapproximate(projection_matrix * Rsets[i], oa_arg)
        end
    end

    return RsetsProj
end
