struct GLGM06 <: ContinuousPost
    options::TwoLayerOptions

    function GLGM06(𝑂::Options)
        normalized_𝑂 = validate_and_wrap_options(𝑂, options_GLGM06();
            validation=validation_GLGM06,
            normalization=normalization_GLGM06!)
        return new(normalized_𝑂)
    end
end

# convenience constructor from pairs of symbols
GLGM06(𝑂::Pair{Symbol,<:Any}...) = GLGM06(Options(Dict{Symbol,Any}(𝑂)))

# default options
GLGM06() = GLGM06(Options())

# out-of-place initialization
init(𝒫::GLGM06, 𝑆::AbstractSystem, 𝑂::Options) = init!(𝒫, 𝑆, copy(𝑂))

# in-place initialization
function init!(𝒫::GLGM06, 𝑆::AbstractSystem, 𝑂::Options)

    # state dimension for (purely continuous or purely discrete systems)
    𝑂[:n] = statedim(𝑆)

    # solver-specific options (adds default values for unspecified options)
    validate_solver_options_and_add_default_values!(𝑂)

    if 𝑂[:project_reachset]
        𝑂[:output_function] = nothing
    else
        𝑂[:output_function] = 𝑂[:projection_matrix]
    end

    return 𝑂
end

function validation_GLGM06(𝑂::TwoLayerOptions)
    return nothing
end

function normalization_GLGM06!(𝑂::TwoLayerOptions)
    return nothing
end

function post(𝒫::GLGM06, 𝑆::AbstractSystem, invariant::LazySet, 𝑂::Options)

    # ==================================
    # Initialization and discretization
    # ==================================
   
    𝑂 = TwoLayerOptions(merge(𝑂, 𝒫.options.specified), 𝒫.options.defaults)
    max_order = 𝑂[:max_order]
    N = round(Int, 𝑂[:T] / 𝑂[:δ])

    𝑆d = discretize(𝑆, 𝑂[:δ], algorithm="forward", set_operations="zonotope")

    Ω0, Φ = 𝑆d.x0, 𝑆d.s.A

    # =====================
    # Flowpipe computation
    # =====================

    # preallocate output
    RSets = Vector{Zonotope}(undef, N)

    info("Reachable States Computation...")
    @timing begin
    if inputdim(𝑆d) == 0 # homogeneous system
        reach_homog!(RSets, Ω0, Φ, N, max_order)
    else
        error("not implemented")
        #=
        # inputs contain the origin
        if zeros(𝑂[:n]) ∈ next_set(𝑈)
            Rsets = reach_inhomog_case1(𝑆, invariant, 𝑂)
        else
            Rsets = reach_inhomog_case2(𝑆, invariant, 𝑂)
        end
        =#
    end
    end # timing

    # ===========
    # Projection
    # ===========
    if 𝑂[:project_reachset] || 𝑂[:projection_matrix] != nothing
        info("Projection...")
        RsetsProj = @timing project(RSets, 𝑂)
    else
        RsetsProj = Rsets
    end

    return ReachSolution(RsetsProj, 𝑂)
end

# ===============================================================
# Homogeneous case
# ===============================================================
function reach_homog!(HR::Vector{Zonotope},
                      Ω0::Zonotope,
                      Φ::AbstractMatrix,
                      N::Int,
                      max_order::Int)
    
    HR[1] = Ω0
    k = 1
    while k < N
        HR_next = linear_map(Φ, HR[k])
        if order(HR_next) > max_order
            HR_next = reduce_order(HR_next, max_order)
        end
        HR[k+1] = HR_next
        k = k + 1
    end
    return HR
end

# ===============================================================
# Inhomogeneous case
# ===============================================================
