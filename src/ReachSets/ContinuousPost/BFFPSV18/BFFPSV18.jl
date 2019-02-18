export BFFPSV18

# ===============================================================
# Bogomolov, Forets, Frehse, Podelski, Schilling, Viry. HSCC 2018
# ===============================================================

# dummy functions for option :lazy_inputs_interval
lazy_inputs_interval_always = (k -> true)
lazy_inputs_interval_never = (k -> false)

function ispartition(partition::AbstractVector{<:AbstractVector{Int}})
    current = 1
    for block in partition
        for i in block
            if i != current
                return false
            end
            current += 1
        end
    end
    return true
end

function options_BFFPSV18()
    return OptionSpec[
        # general options
        OptionSpec(:approx_model, "forward", domain=String, domain_check=(
            v  ->  v in ["forward", "backward", "firstorder", "nobloating"]),
            info="model for bloating/continuous time analysis"),
        OptionSpec(:algorithm, "explicit", domain=String, domain_check=(
            v  ->  v in ["explicit", "wrap"]), info="algorithm backend"),
        OptionSpec(:δ, 1e-2, domain=Float64, aliases=[:sampling_time],
            domain_check=(v  ->  v > 0.), info="time step"),
        OptionSpec(:vars, Int[], domain=AbstractVector{Int}, domain_check=(
            v  ->  length(v) > 0 && all(e -> e > 0, v)),
            info="variables of interest; default: all variables"),
        OptionSpec(:partition, [Int[]],
            domain=AbstractVector{<:AbstractVector{Int}}, domain_check=
            ispartition,
            info="block partition; a block is represented by a vector " *
                 "containing its indices"),

        # discretization options
        OptionSpec(:lazy_sih, false, domain=Bool,
            info="use a lazy symmetric interval hull in discretization?"),
        OptionSpec(:lazy_expm, false, domain=Bool,
            info="use a lazy matrix exponential all the time?"),
        OptionSpec(:lazy_expm_discretize, false, domain=Bool,
            info="use a lazy matrix exponential in discretization?"),
        OptionSpec(:pade_expm, false, domain=Bool,
            info="use the Padé approximant method (instead of Julia's " *
                 " built-in 'exp') to compute the lazy matrix exponential " *
                 "in discretization?"),
        OptionSpec(:assume_sparse, false, domain=Bool,
            info="use an analysis for sparse discretized matrices?"),

        # reachability options
        OptionSpec(:lazy_X0, false, domain=Bool,
            info="keep the discretized and decomposed initial states a lazy " *
                 "set?"),
        OptionSpec(:lazy_inputs_interval, lazy_inputs_interval_always,
            domain=Union{Int, Function},
            domain_check=(v  ->  !(v isa Int) || v >= -1),
            info="length of interval in which the inputs are handled as a " *
                 "lazy set (``-1`` for 'never'); may generally also be a " *
                 "predicate over indices; the default corresponds to ``-1``"),

        # approximation options
        OptionSpec(:block_types, nothing, domain=Union{Nothing,
            Dict{Type{<:LazySet}, AbstractVector{<:AbstractVector{Int}}}},
            info="short hand to set ':block_types_init' and " *
                 "':block_types_iter'"),
        OptionSpec(:block_types_init, nothing, domain=Union{Nothing,
            Dict{Type{<:LazySet}, AbstractVector{<:AbstractVector{Int}}}},
            info="set type for the approximation of the initial states for " *
                 "each block"),
        OptionSpec(:block_types_iter, nothing, domain=Union{Nothing,
            Dict{Type{<:LazySet}, AbstractVector{<:AbstractVector{Int}}}},
            info="set type for the approximation of the states ``X_k``, " *
                 "``k>0``, for each block"),
        OptionSpec(:ε, Inf, domain=Float64, domain_check=(v  ->  v > 0.),
            info="short hand to set `:ε_init` and `:ε_iter`"),
        OptionSpec(:ε_init, Inf, domain=Float64, domain_check=(v  ->  v > 0.),
            info="error bound for the approximation of the initial states" *
                 "(during decomposition)"),
        OptionSpec(:ε_iter, Inf, domain=Float64, domain_check=(v  ->  v > 0.),
            info="error bound for the approximation of the states ``X_k``, " *
                 "``k>0``"),
        OptionSpec(:set_type, Hyperrectangle, domain=Union{Type{HPolygon},
            Type{Hyperrectangle}, Type{LazySets.Interval}},
            info="short hand to set `:set_type_init` and `:set_type_iter`"),
        OptionSpec(:set_type_init, Hyperrectangle, domain=Union{Type{HPolygon},
            Type{Hyperrectangle}, Type{LazySets.Interval}},
            info="set type for the approximation of the initial states" *
                 "(during decomposition)"),
        OptionSpec(:set_type_iter, Hyperrectangle, domain=Union{Type{HPolygon},
            Type{Hyperrectangle}, Type{LazySets.Interval}},
            info="set type for the approximation of the states ``X_k``, " *
                 "``k>0``"),
        OptionSpec(:template_directions, :nothing, domain=Symbol,
            domain_check=(v::Symbol  ->  v in [:box, :oct, :boxdiag, :nothing]),
            info="short hand to set `template_directions_init` and " *
                 "`template_directions_iter`"),
        OptionSpec(:template_directions_init, :nothing, domain=Symbol,
            domain_check=(v::Symbol  ->  v in [:box, :oct, :boxdiag, :nothing]),
            info="directions to use for the approximation of the initial " *
                 "states (during decomposition)"),
        OptionSpec(:template_directions_iter, :nothing, domain=Symbol,
            domain_check=(v::Symbol  ->  v in [:box, :oct, :boxdiag, :nothing]),
            info="directions to use for the approximation of the states " *
                 "``X_k``, ``k>0``, for each block"),

        # convenience options
        OptionSpec(:assume_homogeneous, false, domain=Bool,
            info="ignore dynamic inputs during the analysis?"),
        OptionSpec(:eager_checking, true, domain=Bool,
            info="terminate as soon as property violation was detected?"),
    ]
end

function normalization_BFFPSV18!(𝑂::TwoLayerOptions)
    # :lazy_inputs_interval option: convert integers to functions
    if haskey_specified(𝑂, :lazy_inputs_interval)
        v = 𝑂[:lazy_inputs_interval]
        if v isa Int
            if v == -1
                𝑂.specified[:lazy_inputs_interval] = lazy_inputs_interval_never
            elseif v == 0
                𝑂.specified[:lazy_inputs_interval] = lazy_inputs_interval_always
            else
                𝑂.specified[:lazy_inputs_interval] = (k -> k % v == 0)
            end
        end
    end

    # :block_types options
    block_types = nothing
    dict_type = Dict{Type{<:LazySet}, AbstractVector{<:AbstractVector{Int}}}
    if !haskey_specified(𝑂, :block_types) && haskey(𝑂, :set_type) &&
            haskey_specified(𝑂, :partition)
        𝑂.specified[:block_types] = dict_type(𝑂[:set_type] => copy(𝑂[:partition]))
    end
    if !haskey_specified(𝑂, :block_types_init) && block_types != nothing
        𝑂.specified[:block_types_init] = block_types
    end
    if !haskey_specified(𝑂, :block_types_iter) && block_types != nothing
        𝑂.specified[:block_types_iter] = block_types
    end

    # :ε, :set_type, and :template_directions options
    ε = 𝑂[:ε]
    if haskey_specified(𝑂, :set_type)
        # use the provided set type
        set_type = 𝑂[:set_type]
    elseif ε < Inf
        # use polygons
        set_type = HPolygon
        𝑂[:set_type] = HPolygon
    else
        # use hyperrectangles
        set_type = 𝑂[:set_type]
    end
    #
    if !haskey_specified(𝑂, :ε_init)
        𝑂.specified[:ε_init] =
            (haskey_specified(𝑂, :set_type_init) && 𝑂[:set_type_init] == HPolygon) ||
            (!haskey_specified(𝑂, :set_type_init) && set_type == HPolygon) ?
                ε :
                Inf
    end
    #
    if !haskey_specified(𝑂, :set_type_init)
        𝑂.specified[:set_type_init] = 𝑂[:ε_init] < Inf ? HPolygon : set_type
    end
    #
    if !haskey_specified(𝑂, :template_directions_init)
        𝑂.specified[:template_directions_init] =
            haskey_specified(𝑂, :template_directions_init) ?
                𝑂[:template_directions_init] :
                haskey_specified(𝑂, :template_directions) ?
                    𝑂[:template_directions] :
                    :nothing
    end
    #
    if !haskey_specified(𝑂, :ε_iter)
        𝑂.specified[:ε_iter] =
            (haskey_specified(𝑂, :set_type_iter) && 𝑂[:set_type_iter] == HPolygon) ||
            (!haskey_specified(𝑂, :set_type_iter) && set_type == HPolygon) ?
                ε :
                Inf
    end
    #
    if !haskey_specified(𝑂, :set_type_iter)
        𝑂.specified[:set_type_iter] = 𝑂[:ε_iter] < Inf ? HPolygon : set_type
    end
    #
    if !haskey_specified(𝑂, :template_directions_iter)
        𝑂.specified[:template_directions_iter] =
            haskey_specified(𝑂, :template_directions_iter) ?
                𝑂[:template_directions_iter] :
                haskey_specified(𝑂, :template_directions) ?
                    𝑂[:template_directions] :
                    :nothing
    end
    #

    nothing
end

function validation_BFFPSV18(𝑂)
    # lazy_expm_discretize & lazy_expm
    if !𝑂[:lazy_expm_discretize] && 𝑂[:lazy_expm]
        throw(DomainError(𝑂[:lazy_expm_discretize], "cannot use option " *
            "':lazy_expm' with deactivated option ':lazy_expm_discretize'"))
    end

    # block_types
    if haskey_specified(𝑂, :block_types)
        for (key, value) in 𝑂[:block_types]
            if !(key <: LazySet)
                 throw(DomainError(key, "the keys of the `:block_types` " *
                                        "dictionary should be lazy sets"))
            elseif !(typeof(value) <: AbstractVector{<:AbstractVector{Int}})
                throw(DomainError(value, "the values of the `:block_types` " *
                                         "dictionary should be vectors of " *
                                         "vectors"))
            end
        end
    end

    # ε-close approximation
    if (𝑂[:ε_init] < Inf && 𝑂[:set_type_init] != HPolygon) ||
       (𝑂[:ε_iter] < Inf && 𝑂[:set_type_iter] != HPolygon)
        throw(DomainError("ε-close approximation is only supported with the " *
                          "set type 'HPolygon'"))
    end

    nothing
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
        normalized_𝑂 = validate_and_wrap_options(𝑂, options_BFFPSV18();
            validation=validation_BFFPSV18,
            normalization=normalization_BFFPSV18!)
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

    # :vars option; default: all variables
    if haskey_specified(𝒫.options, :partition)
        𝑂validated[:vars] = 𝒫.options[:vars]
    else
        𝑂validated[:vars] = 1:𝑂validated[:n]
    end

    # :partition option: use 1D blocks
    if haskey_specified(𝒫.options, :partition)
        𝑂validated[:partition] = 𝒫.options[:partition]
    else
        𝑂validated[:partition] = [[i] for i in 1:𝑂validated[:n]]
    end

    # :blocks option (internal only)
    # list of all interesting block indices in the partition
    𝑂validated[:blocks] = compute_blocks(𝑂validated[:vars], 𝑂validated[:partition])

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

function compute_blocks(vars, partition)
    blocks = Vector{Int}()
    sizehint!(blocks, length(vars))
    next = 0
    var_idx = 1
    for (i, block) in enumerate(partition)
        next += length(block)
        if vars[var_idx] <= next
            push!(blocks, i)
            var_idx += 1
            while var_idx <= length(vars) && vars[var_idx] <= next
                var_idx += 1
            end
            if var_idx > length(vars)
                break
            end
        end
    end
    @assert var_idx == length(vars) + 1
    sizehint!(blocks, length(blocks))
    return blocks
end
