"""
Utils module containing helper functions and macros for, e.g., block
decompositions and visualization.
"""
module Utils

using LazySets

# Visualization
export print_sparsity,
       plot_sparsity,
       print_last_interval,
       range_last_x_percent,
       add_plot_labels

# Block Structure
export @block_id,
       add_dimension,
       block_to_set_map,
       convert_partition

# Usability
export @filename_to_png,
       @relpath

# Extension of MathematicalSystems for use inside Reachability.jl
include("systems.jl")

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
    return quote
        f = @__FILE__
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
    block_to_set_map(dict::Dict{Type{<:LazySet},
                                AbstractVector{<:AbstractVector{Int}}})

Invert a map (set type -> block structure) to a map (block index -> set type).

### Input

- `dict` -- map (set type -> block structure)

### Ouput

Vector mapping a block index to the respective set type.
"""
function block_to_set_map(dict::Dict{Type{<:LazySet},
                                     AbstractVector{<:AbstractVector{Int}}})
    # we assume that blocks do not overlap
    set_type = Vector{Type{<:LazySet}}()
    initial_block_indices = Vector{Int}()
    @inbounds for (key, val) in dict
        for bi in val
            push!(set_type, key)
            push!(initial_block_indices, bi[1])
        end
    end
    s = sortperm(initial_block_indices)
    return set_type[s]
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
"""
function convert_partition(partition::AbstractVector{<:AbstractVector{Int}})::Union{
        Vector{Int},
        Vector{UnitRange{Int}},
        Vector{Union{UnitRange{Int}, Int}}}
    # are there 1D blocks and/or kD blocks?
    has_1D_blocks = false
    has_kD_blocks = false
    for (i, block) in enumerate(partition)
        if length(block) == 1
            has_1D_blocks = true
            if has_kD_blocks
                break
            end
        else
            has_kD_blocks = true
            if has_1D_blocks
                break
            end
        end
    end

    # use optimal partition type
    if !has_1D_blocks
        partition_out = Vector{UnitRange{Int}}(length(partition))
    elseif !has_kD_blocks
        partition_out = Vector{Int}(length(partition))
    else
        partition_out = Vector{Union{UnitRange{Int}, Int}}(length(partition))
    end

    # convert old partition representation to new one
    for (i, block) in enumerate(partition)
        if length(block) == 1
            partition_out[i] = block[1]
        else
            partition_out[i] = block[1]:block[end]
        end
    end
    return partition_out
end

end # module
