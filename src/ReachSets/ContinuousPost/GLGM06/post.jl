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
