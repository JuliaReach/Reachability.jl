# ===============================================================
# Homogeneous case
# ===============================================================
function reach_homog!(HR::Vector{ReachSet{Zonotope{Float64}, Float64}},
                      Ω0::Zonotope,
                      Φ::AbstractMatrix,
                      N::Int,
                      δ::Float64,
                      max_order::Int)

    # save timestamps with the reach set
    t0, t1 = zero(δ), δ

    # initial reach set
    HR[1] = ReachSet(Ω0, t0, t1)

    k = 2
    while k <= N
        HR_next = linear_map(Φ, HR[k-1].X)
        HR_next_red = reduce_order(HR_next, max_order)
        t0 = t1; t1 += δ
        HR[k] = ReachSet(HR_next_red, t0, t1)
        k += 1
    end
    return HR
end

# ===============================================================
# Inhomogeneous case
# ===============================================================
function reach_inhomog!(X::Vector{ReachSet{Zonotope{Float64}, Float64}},
                        Ω0::Zonotope,
                        U::ConstantInput,
                        Φ::AbstractMatrix,
                        N::Int,
                        δ::Float64,
                        max_order::Int)

    # initialization
    t0, t1 = zero(δ), δ
    V = next_set(U)

    X[1] = ReachSet(Ω0, t0, t1)
    Wk₊ = V
    Φ_power_k = copy(Φ)
    Φ_power_k_cache = similar(Φ)

    # update
    k = 2
    while k <= N
        Xk = minkowski_sum(linear_map(Φ_power_k, Ω0), Wk₊)
        Wk₊ = minkowski_sum(Wk₊, linear_map(Φ_power_k, V))

        t0 = t1; t1 += δ
        Xk_red = reduce_order(Xk, max_order)
        X[k] = ReachSet(Xk_red, t0, t1)

        Wk₊ = reduce_order(Wk₊, max_order)

        _A_mul_B!(Φ_power_k_cache, Φ_power_k, Φ)
        copyto!(Φ_power_k, Φ_power_k_cache)
        k += 1
    end

    return X 
end
