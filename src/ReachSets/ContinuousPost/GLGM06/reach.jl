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
function reach_inhomog!(X::Vector{ReachSet{Zonotope, Float64}},
                        Ω0::Zonotope,
                        U::ConstantInput,
                        Φ::AbstractMatrix,
                        N::Int,
                        δ::Float64,
                        max_order::Int)

    # preallocation
    W = Vector{Zonotope}(undef, N+1)

    # initialization
    k = 1
    n = size(Φ, 1)
    t0, t1 = zero(δ), δ
    X[k] = ReachSet{Zonotope, Float64}(Ω0, t0, t1)
    W[k] = ZeroSet(n)
    V = next_set(U)
    W[k+1] = V
    Φ_power_k = copy(Φ)

    # update
    k += 1
    while k <= N
        Xk = minkowski_sum(linear_map(Φ_power_k, Ω0), W[k])
        Wk₊ = minkowski_sum(W[k], linear_map(Φ_power_k, V))

        t0 = t1; t1 += δ
        if order(Xk) > max_order
            Xk = reduce_order(Xk, max_order)
        end
        X[k] = ReachSet{Zonotope, Float64}(Xk, t0, t1)

        if order(Wk₊) > max_order
            Wk₊ = reduce_order(Wk₊, max_order)
        end
        W[k+1] = Wk₊

        Φ_power_k = Φ_power_k * Φ
        k += 1
    end

    return X 
end
