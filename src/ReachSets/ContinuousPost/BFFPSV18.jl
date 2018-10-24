export BFFPSV18

# ===============================================================
# Bogomolov, Forets, Frehse, Podelski, Schilling, Viry. HSCC 2018
# ===============================================================

"""
    BFFPSV18 <: ContinuousPost

Implementation of the reachability algorithm for purely continuous linear affine
systems using block decompositons by
Bogomolov, Forets, Frehse, Podelski, Schilling, Viry. HSCC 2018

### Fields

- `options` -- an `Options` structure that holds the algorithm-specific options

### Algorithm

See [Reach Set Approximation through Decomposition with Low-dimensional Sets
and High-dimensional Matrices](https://dl.acm.org/citation.cfm?id=3178128).
S. Bogomolov, M. Forets, G. Frehse, A. Podelski, C. Schilling, F. Viry.
HSCC '18 Proceedings of the 21st International Conference on Hybrid Systems:
Computation and Control (part of CPS Week).
"""
struct BFFPSV18 <: ContinuousPost
    options::Options

    function BFFPSV18(𝑂::Options)
        𝑂copy = copy(𝑂)
        # TODO set defaults for this algorithm here
        return new(𝑂copy)
    end
end

# convenience constructor from pairs of symbols
BFFPSV18(𝑂::Pair{Symbol,<:Any}...) = BFFPSV18(Options(Dict{Symbol,Any}(𝑂)))

# default options
BFFPSV18() = BFFPSV18(Options())

init(𝒫::BFFPSV18, 𝒮::AbstractSystem, 𝑂::Options) = init!(𝒫, 𝒮, copy(𝑂))

function init!(𝒫::BFFPSV18, 𝒮::AbstractSystem, 𝑂::Options)
    # state dimension for (purely continuous or purely discrete systems)
    𝑂copy = copy(𝑂)
    𝑂copy[:n] = statedim(𝒮)

    # solver-specific options (adds default values for unspecified options)
    𝑂validated = validate_solver_options_and_add_default_values!(𝑂copy)

    # Input -> Output variable mapping
    𝑂validated[:inout_map] = inout_map_reach(𝑂validated[:partition], 𝑂validated[:blocks], 𝑂validated[:n])

    if 𝑂validated[:project_reachset]
        𝑂validated[:output_function] = nothing
    else
        𝑂validated[:output_function] = 𝑂validated[:projection_matrix]
    end

    return 𝑂validated
end

"""
    post(𝒫::BFFPSV18, 𝒮::AbstractSystem, invariant, 𝑂::Options)

Calculate the reachable states of the given initial value problem using `BFFPSV18`.

### Input

- `𝒫` -- post operator of type `BFFPSV18`
- `𝒮` -- sytem, initial value problem for a continuous ODE
- `invariant` -- constraint invariant on the mode
- `𝑂` -- algorithm-specific options
"""
function post(𝒫::BFFPSV18, 𝒮::AbstractSystem, invariant, 𝑂::Options)
    # convert matrix
    system = matrix_conversion(𝒮, 𝑂)

    if 𝑂[:mode] == "reach"
        info("Reachable States Computation...")
        tic()
        Rsets = reach(𝒮, invariant, 𝑂)
        info("- Total")
        tocc()

        # Projection
        if 𝑂[:project_reachset] || 𝑂[:projection_matrix] != nothing
            info("Projection...")
            tic()
            RsetsProj = project(Rsets, 𝑂)
            tocc()
        else
            RsetsProj = Rsets
        end

        return ReachSolution(RsetsProj, 𝑂)

    elseif 𝑂[:mode] == "check"
        info("invariants are currently not supported in 'check' mode")

        # Input -> Output variable mapping in property
        𝑂[:property] = inout_map_property(𝑂[:property], 𝑂[:partition], 𝑂[:blocks], 𝑂[:n])

        # =================
        # Property checking
        # =================
        info("Property Checking...")
        tic()
        answer = check_property(𝒮, 𝑂)
        info("- Total")
        tocc()

        if answer == 0
            info("The property is satisfied!")
            return CheckSolution(true, -1, 𝑂)
        else
            info("The property may be violated at index $answer," *
                " (time point $(answer * options[:δ]))!")
            return CheckSolution(false, answer, 𝑂)
        end
    else
        error("unsupported mode $(options[:mode])")
    end # mode
end
