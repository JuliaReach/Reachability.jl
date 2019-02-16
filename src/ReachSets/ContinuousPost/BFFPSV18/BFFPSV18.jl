export BFFPSV18

# ===============================================================
# Bogomolov, Forets, Frehse, Podelski, Schilling, Viry. HSCC 2018
# ===============================================================

function options_BFFPSV18()
    return OptionSpec[
        OptionSpec(:approx_model, "forward", domain=String, domain_check=(
            v  ->  v in ["forward", "backward", "firstorder", "nobloating"]),
            info="model for bloating/continuous time analysis"),
    ]
end

"""
    BFFPSV18 <: ContinuousPost

Implementation of the reachability algorithm for purely continuous linear
time-invariant systems using block decompositons by S. Bogomolov, M. Forets,
G. Frehse, A. Podelski, C. Schilling and F. Viry [1].

### Fields

- `options` -- an `Options` structure that holds the algorithm-specific options

### Notes

The following options are available:

```julia
$(print_option_spec(options_BFFPSV18()))
```

### Algorithm

We refer to [1] for technical details.

[1] [Reach Set Approximation through Decomposition with Low-dimensional Sets
and High-dimensional Matrices](https://dl.acm.org/citation.cfm?id=3178128).
S. Bogomolov, M. Forets, G. Frehse, A. Podelski, C. Schilling, F. Viry.
HSCC '18 Proceedings of the 21st International Conference on Hybrid Systems:
Computation and Control (part of CPS Week).
"""
struct BFFPSV18 <: ContinuousPost
    options::TwoLayerOptions

    function BFFPSV18(𝑂::Options)
        normalized_𝑂 = validate_and_wrap_options(𝑂, options_BFFPSV18())
        return new(normalized_𝑂)
    end
end

# convenience constructor from pairs of symbols
BFFPSV18(𝑂::Pair{Symbol,<:Any}...) = BFFPSV18(Options(Dict{Symbol,Any}(𝑂)))

# default options
BFFPSV18() = BFFPSV18(Options())

init(𝒫::BFFPSV18, 𝑆::AbstractSystem, 𝑂::Options) = init!(𝒫, 𝑆, copy(𝑂))

function init!(𝒫::BFFPSV18, 𝑆::AbstractSystem, 𝑂::Options)
    # state dimension for (purely continuous or purely discrete systems)
    𝑂copy = copy(𝑂)
    𝑂copy[:n] = statedim(𝑆)

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
    post(𝒫::BFFPSV18, 𝑆::AbstractSystem, invariant, 𝑂::Options)

Calculate the reachable states of the given initial value problem using `BFFPSV18`.

### Input

- `𝒫` -- post operator of type `BFFPSV18`
- `𝑆` -- sytem, initial value problem for a continuous ODE
- `invariant` -- constraint invariant on the mode
- `𝑂` -- algorithm-specific options
"""
function post(𝒫::BFFPSV18, 𝑆::AbstractSystem, invariant, 𝑂::Options)
    # TODO temporary hack for refactoring
    𝑂 = TwoLayerOptions(merge(𝑂, 𝒫.options.specified), 𝒫.options.defaults)

    # convert matrix
    system = matrix_conversion(𝑆, 𝑂)

    if 𝑂[:mode] == "reach"
        info("Reachable States Computation...")
        @timing begin
            Rsets = reach(𝑆, invariant, 𝑂)
            info("- Total")
        end

        # Projection
        if 𝑂[:project_reachset] || 𝑂[:projection_matrix] != nothing
            info("Projection...")
            RsetsProj = @timing project(Rsets, 𝑂)
        else
            RsetsProj = Rsets
        end

        return ReachSolution(RsetsProj, 𝑂)

    elseif 𝑂[:mode] == "check"
        info("invariants are currently not supported in 'check' mode")

        # Input -> Output variable mapping in property
        property = inout_map_property(𝑂[:property], 𝑂[:partition], 𝑂[:blocks], 𝑂[:n])

        # =================
        # Property checking
        # =================
        info("Property Checking...")
        @timing begin
            answer = check_property(𝑆, property, 𝑂)
            info("- Total")
        end

        if answer == 0
            info("The property is satisfied!")
            return CheckSolution(true, -1, 𝑂)
        else
            info("The property may be violated at index $answer," *
                " (time point $(answer * 𝑂[:δ]))!")
            return CheckSolution(false, answer, 𝑂)
        end
    else
        error("unsupported mode $(𝑂[:mode])")
    end # mode
end
