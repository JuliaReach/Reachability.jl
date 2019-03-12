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
