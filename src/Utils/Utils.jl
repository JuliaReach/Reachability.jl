"""
Utils module containing helper functions and macros for, e.g., block
decompositions and visualization.
"""
module Utils

using LazySets, MathematicalSystems, HybridSystems

include("../compat.jl")

import Reachability.@timing

# Visualization
export print_sparsity,
       plot_sparsity,
       print_last_interval,
       range_last_x_percent,
       add_plot_labels

# Block Structure
export @block_id,
       add_dimension,
       convert_partition

# Usability
export @filename_to_png,
       @relpath

# internal conversion
export template_direction_symbols,
       matrix_conversion,
       matrix_conversion_lazy_explicit

# internal normalization
export normalize,
       distribute_initial_set

# Extension of MathematicalSystems for use inside Reachability.jl
include("systems.jl")

# normalization
include("normalization.jl")

# abstract solution type
include("AbstractSolution.jl")

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

Prints the max/min values for the last element in the array of two-dimensional
sets in the second dimension.

### Input

- `RsetsProj` -- (projected) reach set given as an array of two-dimensional sets
"""
function print_last_interval(RsetsProj::Vector{<:LazySet})
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
            if iszero(F(bi, bj))
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
            if iszero(block_row_bi[1:2, (2*bj-1):(2*bj)])
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
    add_plot_labels(plot_vars, [project_output], [plot_labels], [use_subindices])

Creates the axis labels for plotting.

### Input

- `plot_vars`      -- variable indices for plotting (0 stands for time)
- `project_output` -- (optional, default: false) flag indicating if the y axis
                      is an output function
- `plot_labels`    -- (optional, default: empty strings) vector of plot labels,
                      empty strings are ignored
- `use_subindices` -- (optional, default: `true`) if `true`, use subindices
                      for the labels, e.g. `x1` is displayed as `x₁`

### Output

Axis labels with the specified properties.
"""
function add_plot_labels(plot_vars::Vector{Int64};
                         project_output::Bool=false,
                         plot_labels::Vector{String}=["", ""],
                         use_subindices::Bool=true)
    labels = copy(plot_labels)

    # convert integer to subindex, such that e.g. x1 is displayed as x₁
    if use_subindices
        tosubindex = i -> join(Char.(0x2080 .+ convert.(UInt16, digits(i)[end:-1:1])))
    else
        tosubindex = i -> string(i)
    end
    if isempty(labels[1])
        xaxis = plot_vars[1]
        labels[1] = (xaxis == 0) ? "t" : "x" * tosubindex(xaxis)
    end
    if isempty(labels[2])
        yaxis = plot_vars[2]
        labels[2] = (yaxis == 0) ? "t" : (project_output ? "y" : "x" * tosubindex(yaxis))
    end
    return labels
end

# Taken from an older Julia version, see Reachability #414 or
# https://github.com/JuliaLang/julia/blob/788ce677f6043b52caf97323f9de4d6975882561/base/loading.jl#L485
function source_path(default::Union{AbstractString, Nothing}="")
    t = current_task()
    while true
        s = t.storage
        if !isa(s, Nothing) && haskey(s, :SOURCE_PATH)
            return s[:SOURCE_PATH]
        end
        if isa(t, t.parent)
            return default
        end
        t = t.parent
    end
end

"""
   @relpath(name)

Returns the path of the current code file.

### Input

- `name` -- file name

### Output

A string.

### Notes

This macro is handy, e.g., when calling a function that uses a file relative to
its own location.

For example, suppose that the folder `/home/projects/P1` contains the script `s.jl`
which loads the data `x.dat` and performs some calculations with it:

```julia
# file /home/projects/P1/s.jl
d = open(@relpath "x.dat")
# do stuff with d
```
In this example, the macro `@relpath "s.jl"` evaluates to the string
`/home/projects/P1/x.dat`. The script `s.jl` can be evaluated from the command line
from any folder, e.g. suppose that the path working directory is `/home/projects`,
then one can do `julia -e 'include("P1/s.jl")'` without having to change `s.jl`.
"""
macro relpath(name::String)
    return quote
        f = source_path()
        if f == nothing
            pathdir = ""
        else
            pathdir = join(split(f, "/")[1:end-1], "/")
        end
        if !isempty(pathdir)
            pathdir = pathdir * "/"
        end
        pathdir * $name
    end
end

"""
    convert_partition(partition::AbstractVector{<:AbstractVector{Int}})::Union{
        Vector{Int}, Vector{UnitRange{Int}}, Vector{Union{UnitRange{Int}, Int}}}

Convert a partition data structure to a more efficient representation where a 1D
block is represented by the index and a higher-dimensional block is represented
by a `UnitRange{Int}`.

### Input

- `partition` -- partition data structure (a vector of integer vectors)

### Output

New partition representation.

### Notes

We cannot mix `Int` and `UnitRange{Int}`, so we only create `Int` partitions if
all blocks are 1D.
"""
function convert_partition(partition::AbstractVector{<:AbstractVector{Int}})::Union{
        Vector{Int},
        Vector{UnitRange{Int}}}
    # are there kD blocks for k > 1?
    has_kD_blocks = false
    for (i, block) in enumerate(partition)
        if length(block) > 1
            has_kD_blocks = true
            break
        end
    end

    # use optimal partition type
    if !has_kD_blocks
        # only 1D blocks
        partition_out = Vector{Int}(undef, length(partition))
        for (i, block) in enumerate(partition)
            partition_out[i] = block[1]
        end
    else
        # at least one kD block for k > 1
        partition_out = Vector{UnitRange{Int}}(undef, length(partition))
        for (i, block) in enumerate(partition)
            partition_out[i] = block[1]:block[end]
        end
    end

    return partition_out
end

template_direction_symbols = Dict(
    :box     => Approximations.BoxDirections,
    :oct     => Approximations.OctDirections,
    :boxdiag => Approximations.BoxDiagDirections
    )

# sparse/dense matrix conversion
function matrix_conversion(Δ, options; A_passed=nothing)
    if A_passed == nothing
        A = Δ.s.A
        create_new_system = false
    else
        A = A_passed
        create_new_system = true
    end

    if options[:assume_sparse]
        if A isa SparseMatrixExp
            # ignore this case
            A_new = A
        elseif !hasmethod(sparse, Tuple{typeof(A)})
            info("`assume_sparse` option cannot be applied to a matrix of " *
                 "type $(typeof(A)) and will be ignored")
            A_new = A
        elseif !(A isa AbstractSparseMatrix)
            # convert to sparse matrix
            A_new = sparse(A)
            create_new_system = true
        else
            # A is already sparse
            A_new = A
        end
    else
        A_new = A
    end
    if create_new_system
        # set new matrix
        if hasmethod(inputset, Tuple{typeof(Δ.s)})
            Δ = IVP(CLCDS(A_new, Matrix(1.0I, size(A_new)), nothing, inputset(Δ)), Δ.x0)
        else
            Δ = IVP(LDS(A_new), Δ.x0)
        end
    end
    return Δ
end

end # module
