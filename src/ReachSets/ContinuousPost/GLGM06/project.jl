# ======================================
# Functionality to project the flowpipe
# ======================================
import LazySets.Approximations: project

# add a "time" variable by taking the cartesian product of the flowpipe ℱ with each time lapse
function add_time(ℱ::Vector{ReachSet{Zonotope{Float64}}})
    ℱ_with_time = Vector{ReachSet{Zonotope{Float64}}}(undef, length(ℱ))
    @inbounds for i in eachindex(ℱ)
        t0, t1 = time_start(ℱ[i]), time_end(ℱ[i])
        radius = (t1 - t0)/2.0
        Xi = set(ℱ[i]) × Zonotope([t0 + radius], hcat(radius))
        Xi = convert(Zonotope, Xi)
        ℱ_with_time[i] = ReachSet{Zonotope{Float64}}(Xi, t0, t1)
    end
    return ℱ_with_time
end

function project(sol::ReachSolution{Zonotope{Float64}})
    n = dim(set(first(sol.flowpipes[1].reachsets))) # state space dimension
    options = copy(sol.options)
    πvars = sol.options[:plot_vars] # variables for plotting
    @assert length(πvars) == 2

    proj_flowpipes = Vector{Flowpipe}(undef, length(sol.flowpipes))
    for (j, flowpipe) in enumerate(sol.flowpipes)
        if 0 ∈ πvars
            # add the time variable to the flowpipe (assuming it's not already
            # part of the model)
            ℱ = add_time(flowpipe.reachsets)
            n += 1
            options[:n] += 1 # TODO : remove when option is removed
            πvars = copy(πvars)
            πvars[first(indexin(0, πvars))] = n # time index is added in the end
        else
            ℱ = flowpipe.reachsets
        end

        πℱ = Vector{ReachSet{Zonotope{Float64}}}(undef, length(flowpipe.reach_sets))
        M = sparse([1, 2], πvars, [1.0, 1.0], 2, n)
        for i in eachindex(ℱ)
            t0, t1 = time_start(ℱ[i]), time_end(ℱ[i])
            πℱ_i = linear_map(M, set(ℱ[i]))
            πℱ_i = Zonotope(πℱ_i.center, πℱ_i.generators)
            πℱ_i = reduce_order(πℱ_i, options[:max_order])
            πℱ[i] = ReachSet{Zonotope{Float64}}(πℱ_i, t0, t1)
        end
        proj_flowpipes[j] = πℱ
    end
    return ReachSolution(proj_flowpipes, options)
end
