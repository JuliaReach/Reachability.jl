"""
    validate_plot_options_and_add_default_values!(polygons, options)

Validates that the given plot options are supported and adds default values for all
unspecified options.

INPUT:

- ``polygons`` -- sequence of polygons
- ``options``  -- an `Options` object, a dictionary of options
                  Supported options:
                  - `:plot_indices` (which indices of the sequence to plot)
                  - `:plot_backend` (which backend should be used)
                  - `:gridlines` (should grid lines be used)
                  - `:plot_name` (name of the plot file)
                  - `:plot_vars` (variables to plot)
                  - `:plot_labels` (axis labels)
                  - `:plot_color` (plot color)
"""
function validate_plot_options_and_add_default_values!(polygons::Vector{HPolygon}, options::Options)
    dict = options.dict

    # use default values for unspecified options
    # TODO<notification> print message to user for each default option, depending on verbosity level
    if !haskey(dict, :plot_indices)
        dict[:plot_indices] = 1:length(polygons)
    end
    if !haskey(dict, :plot_backend)
        dict[:plot_backend] = "pyplot_savefig"
    end
    if !haskey(dict, :gridlines)
        dict[:gridlines] = true
    end
    if !haskey(dict, :plot_name)
        dict[:plot_name] = "plot.png"
    end
    if !haskey(dict, :plot_vars)
        dict[:plot_vars] = [0, 1]
    end
    if !haskey(dict, :plot_labels)
        dict[:plot_labels] = add_plot_labels(dict[:plot_vars])
    end
    if !haskey(dict, :plot_color)
        dict[:plot_color] = "blue"
    end

    # validation
    for kv_pair in dict
        key::Symbol = kv_pair.first

        # define type/domain constraints for each known key
	domain_constraints = (v  ->  true)
	if key == :plot_indices
            expected_type = AbstractVector{Int64}
        elseif key == :plot_backend
            expected_type = String
            domain_constraints = (v::String  ->  v in ["", "pyplot_savefig", "pyplot_inline", "gadfly"])
        elseif key == :gridlines
            expected_type = Bool
        elseif key == :plot_name
            expected_type = String
        elseif key == :plot_labels
            expected_type = Vector{String}
            domain_constraints = (v::Vector{String}  ->  length(v) == 2)
        elseif key == :plot_vars
            expected_type = Vector{Int64}
            domain_constraints = (v::Vector{Int64}  ->  length(v) == 2)
        elseif key == :plot_color
            expected_type = String
        else
            error("Unrecognized option '" * string(key) * "' found.")
        end

        value = kv_pair.second
        # check value type
        if !(value isa expected_type)
            error("Option :" * string(key) * " must be of '" * string(expected_type) * "' type.")
        end
        # check value domain constraints
        if !domain_constraints(value)
            error(string(value) * " is not a valid value for option " * string(key) * ".")
        end
    end

    #return options
end

"""
    plot(polygons, [options])

Plots a sequence of polygons with the given options.

INPUT:

- ``polygons`` -- sequence of polygons
- ``options`` -- options for plotting
"""
@recipe function plot_sol(sol::Rsets2DSolution;
                          seriescolor="blue", label="", grid=true, alpha=0.5)
    #                     plot_indices::AbstractVector{Int64}=1:length(sol.polygons))

    #validate_plot_options_and_add_default_values!(sol.polygons, sol.options)

    seriestype := :shape
#=
    if plot_options[:plot_indices]
         idx = plot_options[:plot_indices]
    else
        idx = 1:length(sol.polygons)
    end
=#
    for i in 1:length(sol.polygons)
        vlist = hcat(vertices_list(sol.polygons[i])...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end