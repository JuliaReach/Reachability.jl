function post(𝒜::GLGM06,
              𝑃::InitialValueProblem{<:AbstractContinuousSystem},
              𝑂::Options)::ReachSolution{Zonotope}

    # ==================================
    # Initialization and discretization
    # ==================================

    𝑂 = merge(𝒜.options.defaults, 𝑂, 𝒜.options.specified)
    max_order = 𝑂[:max_order]
    δ, T = 𝑂[:δ], 𝑂[:T]
    N = round(Int, T / δ)

    # compute and unrwap discretized system
    𝑃_discrete = discretize(𝑃, δ, algorithm=𝑂[:discretization], set_operations="zonotope")
    Ω0, Φ = 𝑃_discrete.x0, 𝑃_discrete.s.A

    # =====================
    # Flowpipe computation
    # =====================

    # preallocate output
    RSets = Vector{ReachSet{Zonotope, Float64}}(undef, N)

    info("Reachable States Computation...")
    @timing begin
    if inputdim(𝑃_discrete) == 0
        reach_homog!(RSets, Ω0, Φ, N, δ, max_order)
    else
        U = inputset(𝑃_discrete)
        reach_inhomog!(RSets, Ω0, U, Φ, N, δ, max_order)
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
