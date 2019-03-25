using TaylorModels: validated_integ
using TaylorSeries: set_variables
using LazySets.Approximations: box_approximation

function post(𝒜::TMJets,
              𝑃::InitialValueProblem{<:Union{BBCS, CBBCS, CBBCCS}, <:LazySet},
              𝑂_global::Options)

    # ==================================
    # Initialization
    # ==================================

    𝑂 = merge(𝒜.options.defaults, 𝑂_global, 𝒜.options.specified)

    # system of ODEs
    f! = 𝑃.s.f
    n = statedim(𝑃)

    # initial time and time horizon
    t0 = 0.0
    T = 𝑂[:T]

    # maximum allowed number of steps
    max_steps = 𝑂[:max_steps]

    # unrap algorithm-specific options
    abs_tol, orderQ, orderT = 𝑂[:abs_tol], 𝑂[:orderQ], 𝑂[:orderT]

    # initial sets
    box_x0 = box_approximation(𝑃.x0)
    q0 = center(box_x0)
    δq0 = IntervalBox(low(box_x0)-q0, high(box_x0)-q0)

    # fix the working variables and maximum order in the global
    # parameters struct (_params_TaylorN_)
    set_variables("x", numvars=length(q0), order=2*orderQ)

    # define the property
    if 𝑂[:mode] == "check"
        property = 𝑂[:property]
    elseif 𝑂[:mode] == "reach"
        if haskey(𝑂, :property)
            property = 𝑂[:property]
        else
            property = (t, x) -> true
        end
    end
    # =====================
    # Flowpipe computation
    # =====================
    invariant = stateset(𝑃.s)

    info("Reachable States Computation...")
    @timing begin
        tTM, xTM = validated_integ(f!, q0, δq0, t0, T, orderQ, orderT, abs_tol,
                                   maxsteps=max_steps, check_property=property, invariant=invariant)
    end

    # convert to hyperrectangle and wrap around the reach solution
    N = length(xTM)
    Rsets = Vector{ReachSet{Hyperrectangle{Float64}, Float64}}(undef, N-1)
    @inbounds for i in 1:N-1
        Hi = convert(Hyperrectangle, xTM[i])
        t0 = tTM[i]; t1 = tTM[i+1]
        Rsets[i] = ReachSet(Hi, t0, t1)
    end

    Rsol = ReachSolution(Rsets, 𝑂)

    # ===========
    # Projection
    # ===========

    if 𝑂[:project_reachset]
        info("Projection...")
        Rsol = @timing project(Rsol)
    end

    return Rsol
end
