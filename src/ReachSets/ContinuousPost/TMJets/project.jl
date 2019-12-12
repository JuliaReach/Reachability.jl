# ======================================
# Functionality to project the flowpipe
# ======================================
import LazySets.Approximations: project
using LazySets.Approximations: overapproximate

# add a "time" variable by taking the cartesian product of the flowpipe ℱ with
# each time lapse, given by an interval
function add_time(ℱ::Vector{ReachSet{Hyperrectangle{T, VC, VR}}}) where {T, VC, VR}
    HT = Hyperrectangle{T, VC, VR}
    ST = CartesianProduct{T, HT, Interval{T, IA.Interval{T}}}
    ℱ_with_time = Vector{ReachSet{ST}}(undef, length(ℱ))
    @inbounds for i in eachindex(ℱ)
        t0, t1 = time_start(ℱ[i]), time_end(ℱ[i])
        Xi = set(ℱ[i]) × Interval(t0, t1)
        ℱ_with_time[i] = ReachSet(Xi, t0, t1)
    end
    return ℱ_with_time
end

function project(sol::ReachSolution{Hyperrectangle{T, VC, VR}}) where {T, VC, VR}
    n = dim(set(first(sol.flowpipes[1].reachsets))) # state space dimension
    options = copy(sol.options)
    πvars = sol.options[:plot_vars] # variables for plotting
    @assert length(πvars) == 2

    proj_flowpipes = Vector{Flowpipe}(undef, length(sol.flowpipes))
    for (j, flowpipe) in enumerate(sol.flowpipes)
        if 0 ∈ πvars # if we want to project along time
            ℱ = add_time(flowpipe.reachsets)
            n += 1
            options[:n] += 1 # TODO : remove when option is removed
            πvars = copy(πvars)
            πvars[first(indexin(0, πvars))] = n # time index is added in the end
        else
            ℱ = flowpipe.reachsets
        end

        HT = Hyperrectangle{T, Vector{T}, Vector{T}} # we hardcode dense vectors here,
        # as using VC and VR would require to adapt overapproximate
        πℱ = Vector{ReachSet{HT}}(undef, length(flowpipe.reach_sets))

        M = sparse([1, 2], πvars, [1.0, 1.0], 2, n)
        for i in eachindex(ℱ)
            t0, t1 = time_start(ℱ[i]), time_end(ℱ[i])
            πℱ_i = overapproximate(M * set(ℱ[i]), Hyperrectangle)
            πℱ[i] = ReachSet(πℱ_i, t0, t1)
        end
        proj_flowpipes[j] = πℱ
    end
    return ReachSolution(proj_flowpipes, options)
end
