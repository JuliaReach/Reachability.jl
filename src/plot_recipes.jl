"""
    plot_sol(sol::ReachSolution; ...)

Plots the solution of a reachability problem in 2D with the given options.

### Input

- `sol`            --  the solution of a reachability problem, projected into
                       two dimensions
- `seriescolor`     -- (optional, default: "blue"): color of the polygons; by
                       default, the same color is used for all of them
- `label`           -- (optional, default: nothing): the legend label
- `grid`            -- (optional, default: true): to use gridlines or not
- `alpha`           -- (optional, default: 0.5): the transparency of the polygons
- `indices`         -- (optional, default: nothing): if given, plot only a sample of
                       the array; this option can be passed through `sol.options`
- `vars`            -- (optional, default: nothing): if given, label the x and y
                       axes with the variables being plotted; this option can be
                       passed through `sol.options`, or as a pair of integers,
                       where 0 indicates the time variable
- `use_subindices`  -- (optional, default: `true`) if `true`, use subindices
                       for the labels, e.g. `x1` is displayed as `x₁`

### Notes

To define your own x and y labels, use the `xlabel` (resp. `ylabel`) keyword
argument. For additional options, consult the Plots.jl reference manual.
"""
@recipe function plot_sol(sol::ReachSolution;
                          seriescolor=:auto,
                          fillcolor=:auto,
                          seriestype=:shape,
                          label="", grid=true, alpha=0.5,
                          indices=nothing, vars=nothing,
                          use_subindices=true)
    @assert dim(sol.Xk[1].X) == 2 "we only support plotting 2D sets"

    options = check_aliases_and_add_default_value(sol)

    if vars != nothing
        vars = add_plot_labels(vars, use_subindices=use_subindices)
        xguide --> vars[1]; yguide --> vars[2]
    elseif options.dict[:plot_vars] != nothing
        vars = add_plot_labels(options.dict[:plot_vars], use_subindices=use_subindices)
        xguide --> vars[1]; yguide --> vars[2]
    end

    if indices == nothing
        indices = options.dict[:plot_indices]
    end

    N = length(indices)
    # rule of thumb for linecolor and fillcolor overlap when plotting many sets
    if N < 300
        linecolor --> :black
    else
        linecolor --> :match
    end

    # Using single list and NaN separators
    vlist = Vector{Vector{Float64}}()
    for i in indices
        append!(vlist, convex_hull(vertices_list(sol.Xk[i].X)))
        push!(vlist, [NaN; NaN])
    end
    vlist = hcat(vlist...)'
    vlist[:, 1], vlist[:, 2]

end

"""
    check_aliases_and_add_default_value(sol::ReachSolution)

Creates a copy of an options structure where aliases have been converted to the
symbol that we use internally.

### Input

- `sol`       -- the solution of a reachability problem, projected in 2D
- `plot_vars` -- (optional) variables for plotting
- `alias`     -- (optional) output_variables
"""
function check_aliases_and_add_default_value(sol::ReachSolution)::Options
    dict = sol.options.dict
    options_copy = Options()
    dict_copy = options_copy.dict

    check_aliases!(dict, dict_copy, [:plot_vars, :output_variables])
    check_aliases_and_add_default_value!(dict, dict_copy, [:plot_indices], 1:length(sol.Xk), false)

    return options_copy
end
