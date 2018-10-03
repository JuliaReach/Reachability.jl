export solve,
       project

function default_algorithm(system::InitialValueProblem)
    algorithm = ""
    s = system.s
    if s isa LinearContinuousSystem ||
       s isa LinearControlContinuousSystem ||
       s isa ConstrainedLinearContinuousSystem ||
       s isa ConstrainedLinearControlContinuousSystem ||
       s isa LinearDiscreteSystem ||
       s isa LinearControlDiscreteSystem ||
       s isa ConstrainedLinearDiscreteSystem ||
       s isa ConstrainedLinearControlDiscreteSystem
       
        algorithm = "BFFPSV18"
    else
        error("no default reachability algorithm available for system of " *
              "type $(typeof(system))")
    end
    return algorithm
end

"""
    solve(system, options)  or  solve(system, :key1 => val1, [...], keyK => valK)

Solves a reachability problem s.t. the given options.
If some options are not defined, we may fall back to default values.

### Input

- `system`    -- a (discrete or continuoues) system specification
- `options`   -- algorithm options for solving the problem
- `algorithm` -- (optional, default: dispatched on the system's type) the
                 reachability algorithm for the computation

### Output

A solution object whose content depends on the input options.

### Notes

To see all available input options, see
`keys(Reachability.available_keywords.dict)`.
"""
function solve(system::InitialValueProblem, options::Options; algorithm::String=default_algorithm(system))
    solve!(system, Options(copy(options.dict)), algorithm=algorithm)
end

solve(system::AbstractSystem, options::Pair{Symbol,<:Any}...) =
    solve(system, Options(Dict{Symbol,Any}(options)))

function solve!(system::InitialValueProblem, options::Options;
                algorithm::String=default_algorithm(system))::AbstractSolution
    if algorithm == "BFFPSV18"
        options = init_BFFPSV18!(system, options)

        # coordinate transformation
        options[:transformation_matrix] = nothing
        if options[:coordinate_transformation] != ""
            info("Transformation...")
            tic()
            (system, transformation_matrix) =
                transform(system, options[:coordinate_transformation])
            tocc()
            options[:transformation_matrix] = transformation_matrix
        end

        # convert matrix
        system = matrix_conversion(system, options)

        if options[:mode] == "reach"
            info("Reachable States Computation...")
            tic()
            Rsets = reach(system, options)
            info("- Total")
            tocc()

            # Projection
            if options[:project_reachset] || options[:projection_matrix] != nothing
                info("Projection...")
                tic()
                RsetsProj = project(Rsets, options)
                tocc()
            else
                RsetsProj = Rsets
            end

            return ReachSolution(RsetsProj, options)

        elseif options[:mode] == "check"

            # Input -> Output variable mapping in property
            options.dict[:property] = inout_map_property(options[:property],
                options[:partition], options[:blocks], options[:n])

            # =================
            # Property checking
            # =================
            info("Property Checking...")
            tic()
            answer = check_property(system, options)
            info("- Total")
            tocc()

            if answer == 0
                info("The property is satisfied!")
                return CheckSolution(true, -1, options)
            else
                info("The property may be violated at index $answer," *
                    " (time point $(answer * options[:δ]))!")
                return CheckSolution(false, answer, options)
            end
        else
            error("unsupported mode $(options[:mode])")
        end # mode
    else
        error("unsupported algorithm $algorithm")
    end # algorithm
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
function project(Rsets::Vector{<:ReachSet}, options::Options)
    plot_vars = copy(options[:plot_vars])
    for i in 1:length(plot_vars)
        if plot_vars[i] != 0
            plot_vars[i] = options[:inout_map][plot_vars[i]]
        end
    end
    reduced_n = sum(x -> x != 0, options[:inout_map])
    output_function = !options[:project_reachset]
    RsetsProj = project_reach(plot_vars,
                              reduced_n,
                              options[:δ],
                              Rsets,
                              options[:algorithm],
                              ε=options[:ε_proj],
                              set_type=options[:set_type_proj],
                              transformation_matrix=options[:transformation_matrix],
                              projection_matrix=options[:projection_matrix],
                              output_function=output_function
                             )
end

project(reach_sol::AbstractSolution) = project(reach_sol.Xk, reach_sol.options)

project(Rsets::Vector{<:ReachSet}, options::Pair{Symbol,<:Any}...) =
    project(Rsets, Options(Dict{Symbol,Any}(options)))

# ===========================================================================
# Bogomolov, Forets, Frehse, Podelski, Schilling, Viry 2018
# ===========================================================================
function init_BFFPSV18!(system, options_input)
    # state dimension for (purely continuous or purely discrete systems)
    options_input.dict[:n] = statedim(system)

    # solver-specific options (adds default values for unspecified options)
    options = validate_solver_options_and_add_default_values!(options_input)

    # Input -> Output variable mapping
    options.dict[:inout_map] = inout_map_reach(options[:partition], options[:blocks], options[:n])

    if options[:project_reachset]
        options[:output_function] = nothing
    else
        options[:output_function] = options[:projection_matrix]
    end

    return options
end
