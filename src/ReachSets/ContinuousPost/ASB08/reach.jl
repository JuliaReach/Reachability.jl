function _add_chunk!(Rsets, Rlin, R̄err)
    @inbounds for i in eachindex(Rlin.Xk)
        Ri = Rlin.Xk[i].X ⊕ R̄err
        Ri = overapproximate(Ri, Hyperrectangle)
        Ri = ReachSet(Ri, Rlin.Xk[i].t_start, Rlin.Xk[i].t_end)
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
    return L̄
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
