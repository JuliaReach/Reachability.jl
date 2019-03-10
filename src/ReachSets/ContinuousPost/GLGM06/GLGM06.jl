export GLGM06

struct GLGM06 <: ContinuousPost
    options::TwoLayerOptions

    function GLGM06(𝑂::Options)
        𝑂new = validate_and_wrap_options(𝑂, options_GLGM06())
        return new(𝑂new)
    end
end

# convenience constructor from pairs of symbols
GLGM06(𝑂::Pair{Symbol,<:Any}...) = GLGM06(Options(Dict{Symbol,Any}(𝑂)))

# default options (they are added in the function validate_and_wrap_options)
GLGM06() = GLGM06(Options())

# out-of-place initialization
init(𝒫::GLGM06, 𝑆::AbstractSystem, 𝑂::Options) = init!(𝒫, 𝑆, copy(𝑂))

function options_GLGM06()

    𝑂spec = Vector{OptionSpec}()

    # step size
    push!(𝑂spec, OptionSpec(:δ, 1e-2, domain=Float64, aliases=[:sampling_time],
                            domain_check=(v  ->  v > 0.), info="time step"))
 
    # discretization
    push!(𝑂spec, OptionSpec(:discretization, "forward", domain=String,
                            info="model for bloating/continuous time analysis"))
            
    push!(𝑂spec, OptionSpec(:sih_method, "concrete", domain=String,
                            info="method to compute the symmetric interval hull in discretization"))

    push!(𝑂spec, OptionSpec(:exp_method, "base", domain=String,
                            info="method to compute the matrix exponential"))

    # approximation options
    push!(𝑂spec, OptionSpec(:max_order, 10, domain=Int,
                            info="maximum allowed order of zonotopes"))

    return 𝑂spec
end

# in-place initialization
function init!(𝒫::GLGM06, 𝑆::AbstractSystem, 𝑂::Options)

    # state dimension
    𝑂[:n] = statedim(𝑆)

    # adds default values for unspecified options
    𝑂init = validate_solver_options_and_add_default_values!(𝑂)

    return 𝑂init
end

function post(𝒫::GLGM06,
              𝑆::AbstractSystem,
              invariant::Union{LazySet, Nothing},
              𝑂::Options)::ReachSolution{Zonotope}

    # ==================================
    # Initialization and discretization
    # ==================================
   
    𝑂 = TwoLayerOptions(merge(𝑂, 𝒫.options.specified), 𝒫.options.defaults)
    max_order = 𝑂[:max_order]
    δ = 𝑂[:δ]
    N = round(Int, 𝑂[:T] / δ)

    # compute and unrwap discretized system
    𝑆d = discretize(𝑆, δ, algorithm=𝑂[:discretization], set_operations="zonotope")
    Ω0, Φ = 𝑆d.x0, 𝑆d.s.A

    # =====================
    # Flowpipe computation
    # =====================

    # preallocate output
    RSets = Vector{ReachSet{Zonotope, Float64}}(undef, N)

    info("Reachable States Computation...")
    @timing begin
    if inputdim(𝑆d) == 0
        reach_homog!(RSets, Ω0, Φ, N, δ, max_order)
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
        RsetsProj = RSets
    end

    return ReachSolution(RsetsProj, 𝑂)
end

# ===============================================================
# Homogeneous case
# ===============================================================
function reach_homog!(HR::Vector{ReachSet{Zonotope, Float64}},
                      Ω0::Zonotope,
                      Φ::AbstractMatrix,
                      N::Int,
                      δ::Float64,
                      max_order::Int)

    # save timestamps with the reach set
    t0, t1 = zero(δ), δ

    # initial reach set
    HR[1] = ReachSet{Zonotope, Float64}(Ω0, t0, t1)

    k = 1
    while k < N
        HR_next = linear_map(Φ, HR[k].X)
        if order(HR_next) > max_order
            HR_next = reduce_order(HR_next, max_order)
        end
        t0 = t1; t1 += δ
        HR[k+1] = ReachSet{Zonotope, Float64}(HR_next, t0, t1)
        k = k + 1
    end
    return HR
end

# ===============================================================
# Inhomogeneous case
# ===============================================================
