"""
Utils module containing helper functions and macros for, e.g., block
decompositions and visualization.
"""
module Utils

using LazySets, ..Systems

# Visualization
export print_sparsity,
       plot_sparsity,
       print_last_interval,
       range_last_x_percent,
       add_plot_labels

# Block Structure
export @block_id,
       add_dimension

# Usability
export @filename_to_png,
       @relpath

"""
    @filename_to_png

Expands into the current filename and transforms the file extension to `png`.

### Examples

If this macro is executed from a script named my_model.jl, then:

```julia
julia> plot_name = @filename_to_png
my_model.png
```
"""
macro filename_to_png()
    return :(split(split(@__FILE__, "/")[end], ".")[1] * ".png")
end

"""
    @block_id v

Returns the block number associated to a given variable/dimension.

It is assumed that the block size is two and that blocks are ordered from top
(first two rows) to bottom (last two rows) in a matrix.

### Examples

```julia
julia> @block_id 4
2
```
"""
macro block_id(v::Int64)
    return :(div($v, 2) + mod($v, 2))
end

"""
    add_dimension(A)

Adds an extra zero row and column to a matrix.

### Input

- `A` -- matrix

### Examples

```julia
julia> A = [0.4 0.25; 0.46 -0.67]
2×2 Array{Float64,2}:
 0.4    0.25
 0.46  -0.67
julia> add_dimension(A)
3×3 Array{Float64,2}:
 0.4    0.25  0.0
 0.46  -0.67  0.0
 0.0    0.0   0.0
```
"""
function add_dimension(A::AbstractMatrix{Float64})::AbstractMatrix{Float64}
    n = size(A, 1)
    return vcat(hcat(A, zeros(n)), zeros(n+1).')
end

"""
    add_dimension(X)

Adds an extra dimension to a LazySet, usually through a Cartesian product.

### Input

- `X` -- a lazy set

### Examples

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
    return X * ZeroSet()
end

function add_dimension(X::ZeroSet)::ZeroSet
    return ZeroSet(dim(X)+1)
end

"""
    add_dimension(cont_sys)

Adds an extra dimension to a continuous system.

### Input

- `cs` -- continuous system

### Examples

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
julia> dim(next_set(sext.U))
4
```
"""
function add_dimension(cs::ContinuousSystem)::ContinuousSystem
    Aext = add_dimension(cs.A)
    X0ext = add_dimension(cs.X0)
    if cs.U isa ConstantNonDeterministicInput
        Uext = add_dimension(next_set(cs.U))
    elseif cs.U isa TimeVaryingNonDeterministicInput
        Uext = Vector{LazySet}(length(cs.U))
        i = 0
        for u in cs.U
            Uext[i+=1] = add_dimension(u)
        end
    else
        error("Unsupported inputs type $(typeof(cs.U)).")
    end
    return ContinuousSystem(Aext, X0ext, Uext)
end

"""
    range_last_x_percent(N, first=0, step=1)

Returns a StepRange for the last x percent of the range 1:N.

### Input

- `N`     -- range length
- `first` -- relative starting index (in percent); a number from [0, 99]
                (default: 0)
- `step`  -- step size (default: 1)

EXAMPLE:

```julia
julia> range_last_x_percent(200, 20, 5)
160:5:200
```
"""
function range_last_x_percent(N::Int64, first::Int64=0, step::Int64=1)
    if first == 0
        start = 1
    elseif first > 0 && first < 100
        start = (100-first)*(floor(Int, N/100))
        start = max(1, start + ((N-start) % step)) # last index should be N
    else
        throw(DomainError("Start index must be in [0, 99]."))
    end
    return start : step : N
end

"""
    print_last_interval(RsetsProj)

Prints the max/min values for the last element in the array of HPolygons in the
second dimension.

### Input

- `RsetsProj` -- (projected) reach set given as an array of HPolygons
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

### Input

- `ϕ`    -- matrix
- `name` -- (optional, default: "") matrix name (for the output message)
"""
function print_sparsity(ϕ::AbstractMatrix{Float64}, name::String="")
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

    print_sparsity_message(name, zero_blocks, all_blocks)
end

function print_sparsity(ϕ::SparseMatrixExp{Float64}, name::String="")
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

    print_sparsity_message(name, zero_blocks, all_blocks)
end

@inline function print_sparsity_message(name::String, zero_blocks::Int, all_blocks::Int)
    sparsity = zero_blocks / all_blocks * 100.
    @printf "sparsity %s: %.3f %% (%d/%d zero blocks)\n" name sparsity zero_blocks all_blocks
end

"""
    plot_sparsity(ϕ, name, [plot_backend])

Plots the sparsity of a matrix ϕ.

### Input

- `ϕ`            -- dense or sparse matrix
- `name`         -- the name of the matrix (for the file name)
- `plot_backend` -- (optional, default: `''`): name of the plotting backend;
                      valid values are:

                 -  `'pyplot_savefig'` -- use PyPlot package, save to a file

                 -  `'pyplot_inline'` -- use PyPlot package, showing in external program

                 - `''` -- (empty string), no plotting
"""
function plot_sparsity(ϕ::Union{AbstractMatrix{Float64}, SparseMatrixExp{Float64}}, name::String, plot_backend::String="")
    if isempty(plot_backend)
        # plot nothing
    elseif plot_backend in ["pyplot_inline", "pyplot_savefig"]
        if !isdefined(:PyPlot)
              error("this backend requires that your script loads the PyPlot module")
        end
        eval(Expr(:using, :PyPlot))
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
        warn("plotting backend $plot_backend is not supported for sparsity plots")
    end
end

"""
    add_plot_labels(plot_vars, [project_output], [plot_labels])

Creates the axis labels for plotting.

### Input

- `plot_vars`      -- variable indices for plotting (0 stands for time)
- `project_output` -- (optional, default: false) flag indicating if the y axis
                      is an output function
- `plot_labels`    -- (optional) vector of plot labels (empty strings are ignored)
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


"""
   @relpath(name)

Returns the path of the current code file.
This is handy, e.g., when calling a function from anywhere that wants to open a
file relative to its own location.

### Input

- `name` -- file name
"""
macro relpath(name::String)
    return :(join(split(@__FILE__, "/")[1:end-1], "/") * "/" * $name)
end

end # module
