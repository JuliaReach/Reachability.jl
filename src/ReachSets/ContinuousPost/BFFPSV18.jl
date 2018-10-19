# ===============================================================
# Bogomolov, Forets, Frehse, Podelski, Schilling, Viry. HSCC 2018
# ===============================================================

struct BFFPSV18 <: ContinuousPost
end

function init(op::BFFPSV18, system, options_input)
    # state dimension for (purely continuous or purely discrete systems)
    options_copy = copy(options_input)
    options_copy.dict[:n] = statedim(system)

    # solver-specific options (adds default values for unspecified options)
    options = validate_solver_options_and_add_default_values!(options_copy)

    # Input -> Output variable mapping
    options.dict[:inout_map] =
        inout_map_reach(options[:partition], options[:blocks], options[:n])

    if options[:project_reachset]
        options[:output_function] = nothing
    else
        options[:output_function] = options[:projection_matrix]
    end

    return options
end

function post(op::BFFPSV18, system, invariant, options)
    # convert matrix
    system = matrix_conversion(system, options)

    if options[:mode] == "reach"
        info("Reachable States Computation...")
        tic()
        Rsets = reach(system, invariant, options)
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
        info("invariants are currently not supported in 'check' mode")

        # Input -> Output variable mapping in property
        options.dict[:property] =
            inout_map_property(options[:property],
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
end
