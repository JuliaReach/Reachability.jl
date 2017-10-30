"""
    plot_sol(sol::Rsets2DSolution; ...)

Plots the solution of a reachability problem in 2D with the given options.

### Input

- `Rsets2DSolution` -- the solution of a reachability problem, projected into
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

### Notes

To define your own x and y labels, use the `xlabel` (resp. `ylabel`) keyword
argument. For additional options, consult the Plots.jl reference manual.
"""
@recipe function plot_sol(sol::Rsets2DSolution;
                          seriescolor="blue", label="", grid=true, alpha=0.5,
                          indices=nothing, vars=nothing)

    seriestype := :shape

    if vars != nothing
        vars = add_plot_labels(vars)
        xguide --> vars[1]; yguide --> vars[2]
    elseif haskey(sol.options.dict, :plot_vars)
        vars = add_plot_labels(sol.options.dict[:plot_vars])
        xguide --> vars[1]; yguide --> vars[2]
    end

    if indices == nothing
        if haskey(sol.options.dict, :plot_indices)
            indices = sol.options.dict[:plot_indices]
        else
            indices = 1:length(sol.polygons)
        end
    end

    for i in indices
        vlist = hcat(vertices_list(sol.polygons[i])...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end
