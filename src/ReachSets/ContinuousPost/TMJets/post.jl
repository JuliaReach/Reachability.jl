using TaylorModels: validated_integ

function post(𝒫::TMJets,
              𝑆::AbstractSystem, # {<:ImplicitContinuousSystem}
              𝑂::Options)::ReachSolution{Zonotope}

    # ==================================
    # Initialization
    # ==================================

    𝑂 = TwoLayerOptions(merge(𝑂, 𝒫.options.specified), 𝒫.options.defaults)

    # system of ODEs
    f! = 𝑆.s
    n = 𝑂[:n]

    # initial and final times, and maximum allowed number of steps
    t0 = 0.0
    T = 𝑂[:T]
    max_steps = 𝑂[:max_steps]

    # unrap algorithm-specific options
    abs_tol, orderQ, orderT = 𝑂[:abs_tol], 𝑂[:orderQ] 𝑂[:orderT]

    # initial sets
    X0 = convert(IntervalBox, 𝑆.x0)
    q0 = mid(X0)
    δq0 = sup.(X0) - mid(X0)

    # returns a TaylorN vector, each entry corresponding to an indep variable
    set_variables("x", numvars=length(q0), order=2*orderQ)

    # define the property
    property = haskey(𝑂, :property) ? 𝑂[:property] : (t, x) -> true

    # =====================
    # Flowpipe computation
    # =====================

    # preallocate output
    RSets = Vector{ReachSet{Hyperrectangle, Float64}}(undef, N)

    info("Reachable States Computation...")
    @timing begin
        tTM, xTM = validated_integ(f!, q0, δq0, t0, T, orderQ, orderT, abs_tol,
                     maxsteps=max_steps, check_property=property)
    end

    # convert to hyperrectangle to wrap around the reach solution
    N = length(xTM)
    RSets = Vector{Hyperrectangle}(undef, N)
    @inbounds for i in eachindex(xTM)
        RSets[i] = convert(Hyperrectangle, xTM[i])
    end

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
