# =======================================
# Functionality to project the flowpipe
# =======================================
import LazySets.Approximations: project

function delete_zero_cols(A::AbstractMatrix)
    nonzero_cols = Vector{Int}()
    for (i, ci) in enumerate(eachcol(A))
        if !iszero(ci)
            push!(nonzero_cols, i)
        end
    end
    return A[:, nonzero_cols]
end

# add a "time" variable by taking the cartesian product of the flowpipe ℱ with each time lapse
function add_time(ℱ)
    ℱ_with_time = Vector{ReachSet{Zonotope{Float64}, Float64}}(undef, length(ℱ))
    @inbounds for i in eachindex(ℱ)
        t0, t1 = ℱ[i].t_start, ℱ[i].t_end
        radius = (t1 - t0)/2.0
        Xi = ℱ[i].X × Zonotope([t0 + radius], hcat(radius))
        Xi = convert(Zonotope, Xi)
        ℱ_with_time[i] = ReachSet{Zonotope{Float64}, Float64}(Xi, t0, t1)
    end
    return ℱ_with_time
end

function project(sol::ReachSolution{Zonotope})
    N = length(sol.Xk)  # number of reach sets
    n = dim(first(sol.Xk).X) # state space dimension
    options = copy(sol.options)
    πℱ = Vector{ReachSet{Zonotope{Float64}, Float64}}(undef, N) # preallocated projected reachsets
    πvars = sol.options[:plot_vars] # variables for plotting
    @assert length(πvars) == 2

    if 0 ∈ πvars
        # add the time variable to the flowpipe (assuming it's not already
        # part of the model)
        ℱ = add_time(sol.Xk)
        n += 1
        options[:n] += 1 # TODO : remove when option is removed
        πvars = copy(πvars)
        πvars[first(indexin(0, πvars))] = n # time index is added in the end
    else
        ℱ = sol.Xk
    end

    M = sparse([1, 2], πvars, [1.0, 1.0], 2, n)
    for i in eachindex(ℱ)
        t0, t1 = ℱ[i].t_start, ℱ[i].t_end
        πℱ_i = linear_map(M, ℱ[i].X)
        πℱ_i = Zonotope(πℱ_i.center, delete_zero_cols(πℱ_i.generators))
        πℱ_i = reduce_order(πℱ_i, options[:max_order])
        πℱ[i] = ReachSet{Zonotope{Float64}, Float64}(πℱ_i, t0, t1)
    end
    return ReachSolution(πℱ, options)
end
