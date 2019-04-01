function linearize(𝑃::IVP{<:BBCS}, δ)

    # nonlinear ODE
    f = 𝑃.s.f
    
    # initial set of current chunk
    𝑋₀i = 𝑃.x0

    # linearization point for current chunk
    c = center(𝑋₀i)
    x̃ = c + δ/2 * f(c)

    # evaluate Jacobian at the linearization point
    Ax̃ = jacobian(f, x̃) #  map(x -> evaluate(x, x̃), Jf)
    bx̃ = f(x̃) - Ax̃ * x̃

    # instantiate linearized system; it doesn't have state constraints
    𝑆lin = ConstrainedAffineContinuousSystem(Ax̃, bx̃, Universe(2));
    𝑃lin = IVP(𝑆lin, 𝑋₀i)
    return x̃, 𝑃lin
end

function _add_chunk!(Rsets::Vector{ReachSet{Hyperrectangle{Float64}, Float64}}, Rlin, R̄err, t0)
    @inbounds for i in eachindex(Rlin.Xk)
        Ri = Rlin.Xk[i].X ⊕ R̄err
        Ri = overapproximate(Ri, Hyperrectangle)
        Ri = ReachSet(Ri, t0 + Rlin.Xk[i].t_start, t0 + Rlin.Xk[i].t_end)
        push!(Rsets, Ri)
    end
    return Rsets
end

function _add_chunk!(Rsets::Vector{ReachSet{Zonotope{Float64}, Float64}}, Rlin, R̄err, t0)
    @inbounds for i in eachindex(Rlin.Xk)
        Ri = minkowski_sum(Rlin.Xk[i].X, convert(Zonotope, R̄err))
        Ri = ReachSet(Ri, t0 + Rlin.Xk[i].t_start, t0 + Rlin.Xk[i].t_end)
        push!(Rsets, Ri)
    end
    return Rsets
end

function admissible_error(Ax̃, δ, θ; n=2)
    @assert n == 2
    Φ₁_Aδ = Φ₁(Ax̃, δ)
    R̄err = Hyperrectangle(zeros(2), θ*δ)
    l̄ = abs.(inv(Φ₁_Aδ)) * θ * δ
    L̄ = Hyperrectangle(zeros(2), l̄)
    return R̄err, L̄
end

function lagrange_remainder(f, Rlin, R̄err, x̃; n=2)
    @assert n == 2
    
    Hf₁ = [∂(f[1], (2, 0)) ∂(f[1], (1, 1));
           ∂(f[1], (1, 1)) ∂(f[1], (0, 2))]
    Hf₂ = [∂(f[2], (2, 0)) ∂(f[2], (1, 1));
           ∂(f[2], (1, 1)) ∂(f[2], (0, 2))]

    R̂lin = ConvexHullArray([Ri.X for Ri in Rlin.Xk]) ⊕ R̄err
    R̂lin_rect = overapproximate(R̂lin, Hyperrectangle)

    ξ = CH(Singleton(x̃), R̂lin_rect)
    ξ_rect = overapproximate(ξ, Hyperrectangle)
    ξ_box = convert(IntervalBox, ξ_rect)

    Hf₁_box = map(Hf_ij -> evaluate(Hf_ij, ξ_box), Hf₁)
    Hf₂_box = map(Hf_ij -> evaluate(Hf_ij, ξ_box), Hf₂)

    R̂lin_zono = convert(Zonotope, R̂lin_rect)

    γ = abs.(R̂lin_zono.center - x̃) + sum(abs.(R̂lin_zono.generators), dims=2)
    
    G = [sup.(abs.(Hf₁_box)), sup.(abs.(Hf₂_box))];
    li_zono = [(1/2 * transpose(γ) * G[i] * γ)[1, 1] for i in 1:n]
    L = Hyperrectangle(zeros(n), li_zono)
    return L
end

