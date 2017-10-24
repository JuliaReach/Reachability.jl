"""
Utils module containing helper functions and macros for block decompositions and
for visualization.
"""
module Utils

using PyPlot, LazySets, ..Systems

# Visualization
export @filename_to_png,
       print_sparsity,
       plot_sparsity,
       print_last_interval,
       range_last_x_percent,
       add_plot_labels

# Block Structure
export @block_id,
       add_dimension

"""
    @filename_to_png

This macro expands into the current filename, and transforms the suffix
to png format.

EXAMPLES:

If this macro is executed from a script name my_model.jl, then:

    julia> plot_name = @filename_to_png
    my_model.png
"""
macro filename_to_png()
    return :(split(split(@__FILE__, "/")[end], ".")[1] * ".png")
end

"""
    @block_id v

Return the block number associated to a given variable.

It is assumed that block size is two and that they are ordered from top (first
row) to bottom (last row).

EXAMPLES:

    julia> @block_id 4
    2
"""
macro block_id(v::Int64)
    return :(div($v, 2) + mod($v, 2))
end

"""
    add_dimension(A)

Add an extra, empty, row and column to a given matrix.

INPUT:

- `A` -- matrix, which can be either dense or sparse

EXAMPLES:

    julia> A = [0.4 0.25; 0.46 -0.67]
    2×2 Array{Float64,2}:
     0.4    0.25
     0.46  -0.67
    julia> add_dimension(A)
    3×3 Array{Float64,2}:
     0.4    0.25  0.0
     0.46  -0.67  0.0
     0.0    0.0   0.0
"""
function add_dimension(A::AbstractArray)::AbstractArray
    n = size(A, 1)
    return vcat(hcat(A, zeros(n)), zeros(n+1).')
end

"""
    add_dimension(X)

Add an extra dimension to a lazy set through a Cartesian product.

INPUT:

- `X` -- a lazy set

EXAMPLES:

```julia
julia> using LazySets
julia> X = BallInf(ones(9), 0.5);
julia> dim(X)
9
julia> Xext = add_dimension(X);
julia> dim(Xext)
10
```
"""
function add_dimension(X::LazySet)::LazySet
    if X isa VoidSet
        return VoidSet(dim(X)+1)
    else
        return X * Singleton([0])
    end
end

"""
    add_dimension(cont_sys)

Add an extra dimension to a continuous system.

INPUT:

- `A` -- matrix, which can be either dense or sparse

EXAMPLES:

```julia
julia> A = sprandn(3, 3, 0.5);
julia> X0 = BallInf(zeros(3), 1.0);
julia> s = ContinuousSystem(A, X0);
julia> sext = add_dimension(s);
julia> Systems.dim(sext)
4
```

If there is an input set, it is also extended:

```julia
julia> U = ConstantNonDeterministicInput(Ball2(ones(3), 0.1));
julia> s = ContinuousSystem(A, X0, U);
julia> sext = add_dimension(s);
julia> dim(start(sext.U).sf)
4
```
"""
function add_dimension(cs::ContinuousSystem)::ContinuousSystem
    Aext = add_dimension(cs.A)
    X0ext = add_dimension(cs.X0)
    if cs.U isa ConstantNonDeterministicInput
        U = start(cs.U).sf
        Uext = add_dimension(U)
        return ContinuousSystem(Aext, X0ext, Uext)
    elseif cs.U isa TimeVaryingNonDeterministicInput
        Uext = LazySet[]
        Ui = start(cs.U)
        push!(Uext, add_dimension(Ui.sf))
        for i in 2:length(cs.U)
            Ui = next(cs.U, Ui)
            push!(Uext, add_dimension(Ui.sf))
        end
        return ContinuousSystem(Aext, X0ext, Uext)
    end
end

"""
    range_last_x_percent(N, first=0, step=1)

Returns a StepRange for the last x percent of the range 1:N.

The user can control the percentage of the start index and the step size.

INPUT:

- ``N``     -- range length
- ``first`` -- relative starting index (in percent); a number from [0, 100)
                (default: 0)
- ``step``  -- step size (default: 1)

EXAMPLE:

    julia> range_last_x_percent(200, 20, 5)
    160:5:200
"""
function range_last_x_percent(N::Int64, first::Int64=0, step::Int64=1)
    if first == 0
        start = 1
    elseif first > 0 && first < 100
        start = (100-first)*(Int)(floor(N/100))
        start = max(1, start + ((N-start) % step)) # last index should be N
    else
        throw(DomainError())
    end
    return start : step : N
end

"""
    print_last_interval(RsetsProj)

Prints the max/min values for the last element in the array of HPolygons.

INPUT:

- ``RsetsProj`` -- (projected) reach set given as an array of HPolygons
"""
function print_last_interval(RsetsProj::Vector{HPolygon})
    last_polygon = RsetsProj[end]
    min_val = -ρ([0., -1.], last_polygon)
    max_val = ρ([0., 1.], last_polygon)

    println("max/min values for the last time step:")
    @printf "[ %.6e ; %.6e ]\n" max_val min_val
end

"""
    print_sparsity(ϕ, name)

Prints the sparsity of a matrix ϕ.

INPUT:

- ``ϕ``    -- a matrix, which can be either dense or sparse
- ``name`` -- the name of the matrix (for the output message)
"""
function print_sparsity(ϕ::Union{Matrix{Float64}, SparseMatrixCSC{Float64, Int64}}, name::String="")
    zero_blocks = 0
    b = (Int64)(size(ϕ, 1) / 2)
    @inline F(bi::Int64, bj::Int64) = ϕ[(2*bi-1):(2*bi), (2*bj-1):(2*bj)]
    for bj in 1:b
        for bi in 1:b
            if findfirst(F(bi, bj)) == 0
                zero_blocks += 1
            end
        end
    end
    all_blocks = b^2
    sparsity = zero_blocks / all_blocks * 100.

    @printf "sparsity %s: %.3f %% (%d/%d zero blocks)\n" name sparsity zero_blocks all_blocks
end

function print_sparsity(ϕ::SparseMatrixExp, name::String="")
    zero_blocks = 0
    b = (Int64)(size(ϕ, 1) / 2)
    @inline F(bi::Int64) = get_rows(ϕ, (2*bi-1):2*bi)
    for bi in 1:b
        block_row_bi = F(bi)
        for bj in 1:b
            if findfirst(block_row_bi[1:2, (2*bj-1):(2*bj)]) == 0
                zero_blocks += 1
            end
        end
    end
    all_blocks = b^2
    sparsity = zero_blocks / all_blocks * 100.

    @printf "sparsity %s: %.3f %% (%d/%d zero blocks)\n" name sparsity zero_blocks all_blocks
end

"""
    plot_sparsity(ϕ, name, [plot_backend])

Plots the sparsity of a matrix ϕ.

INPUT:

- ``ϕ``            -- a matrix, which can be either dense or sparse
- ``name``         -- the name of the matrix (for the file name)
- ``plot_backend`` -- (optional, default: ``''``): name of the plotting backend;
                      valid values are:

                 -  ``'pyplot_savefig'`` -- use PyPlot package, save to a file

                 -  ``'pyplot_inline'`` -- use PyPlot package, showing in external program

                 - ``'gadfly'`` -- use Gadfly package, showing in browser

                 - ``''`` -- (empty string), no plotting
"""
function plot_sparsity(ϕ::Union{Matrix{Float64}, SparseMatrixCSC{Float64, Int64}, SparseMatrixExp}, name::String, plot_backend::String="")
    if isempty(plot_backend)
        # plot nothing
    elseif plot_backend in ["pyplot_inline", "pyplot_savefig"]
        if plot_backend == "pyplot_savefig"
            PyPlot.ioff()  # turn off interactive plotting
            fig = PyPlot.figure()
        end
#       PyPlot.gray()
        if ϕ isa SparseMatrixExp
            explicit_matrix = get_rows(ϕ, 1:size(ϕ, 1))
            PyPlot.spy(explicit_matrix)
        else
            PyPlot.spy(ϕ)
        end
        if plot_backend == "pyplot_savefig"
            PyPlot.savefig(name * "_sparsity.png")
            PyPlot.close(fig)
        end
    else
        @printf "plotting backend %s is not supported for sparsity plots" plot_backend
    end
end

"""
    add_plot_labels(plot_vars, [project_output], [plot_labels])

Creates the axis labels for plotting.

INPUT:

- ``plot_vars``      -- variable indices for plotting (0 stands for time)
- ``project_output`` -- flag indicating if the y axis is an output function
- ``plot_labels``    -- (optional) vector of plot labels (empty strings are ignored)
"""
function add_plot_labels(plot_vars::Vector{Int64}, project_output::Bool=false, plot_labels::Vector{String}=["", ""])
    labels = copy(plot_labels)
    if isempty(labels[1])
        xaxis = plot_vars[1]
        labels[1] = (xaxis == 0) ? "t" : "x" * dec(xaxis)
    end
    if isempty(labels[2])
        yaxis = plot_vars[2]
        labels[2] = (yaxis == 0) ? "t" : (project_output ? "y" : "x" * dec(yaxis))
    end
    return labels
end

end # module
