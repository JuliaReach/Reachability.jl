function post(𝒜::GLGM06,
              𝑃::InitialValueProblem{<:AbstractContinuousSystem},
              𝑂::Options)

    # ==================================
    # Initialization and discretization
    # ==================================

    𝑂 = merge(𝒜.options.defaults, 𝑂, 𝒜.options.specified)
    max_order = 𝑂[:max_order]
    δ, T = 𝑂[:δ], 𝑂[:T]
    N = round(Int, T / δ)

    # compute and unrwap discretized system
    𝑃_discrete = discretize(𝑃, δ, algorithm=𝑂[:discretization],
                                  sih_method=𝑂[:sih_method],
                                  exp_method=𝑂[:exp_method],
                                  set_operations="zonotope")
    Ω0, Φ = 𝑃_discrete.x0, 𝑃_discrete.s.A

    # =====================
    # Flowpipe computation
    # =====================

    # preallocate output
    Rsets = Vector{ReachSet{<:Zonotope{Float64}}}(undef, N)

    info("Reachable States Computation...")
    @timing begin
    if inputdim(𝑃_discrete) == 0
        reach_homog!(Rsets, Ω0, Φ, N, δ, max_order)

    else
        U = inputset(𝑃_discrete)
        reach_inhomog!(Rsets, Ω0, U, Φ, N, δ, max_order)
    end
    end # timing

    Rsol = ReachSolution(Flowpipe(Rsets), 𝑂)

    # ===========
    # Projection
    # ===========

    if 𝑂[:project_reachset]
        info("Projection...")
        Rsol = @timing project(Rsol)
    end

    return Rsol
end
