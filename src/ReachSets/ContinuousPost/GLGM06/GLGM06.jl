export GLGM06

struct GLGM06 <: ContinuousPost
    options::TwoLayerOptions

    function GLGM06(𝑂::Options)
        𝑂new = validate_and_wrap_options(𝑂, options_GLGM06())
        return new(𝑂new)
    end
end

# convenience constructor from pairs of symbols
GLGM06(𝑂::Pair{Symbol,<:Any}...) = GLGM06(Options(Dict{Symbol,Any}(𝑂)))

# default options (they are added in the function validate_and_wrap_options)
GLGM06() = GLGM06(Options())

include("init.jl")
include("post.jl")
include("reach.jl")
include("check.jl")

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
    #return @view(A[:, nonzero_cols])
    return A[:, nonzero_cols]
end

# add a "time" variable by taking the cartesian product of a flowpipe with
# each time lapse
function add_time(sol::ReachSolution{Zonotope})
    N = length(sol.Xk)
    sol_with_time = Vector{ReachSet{Zonotope{Float64}, Float64}}(undef, N)
    @inbounds for i in eachindex(sol.Xk)
        t0, t1 = sol.Xk[i].t_start, sol.Xk[i].t_end
        radius = (t1 - t0)/2.0
        Xk_i = sol.Xk[i].X × Zonotope([t0 + radius], hcat(radius)) # BallInf([t0 + radius], radius)
        Xk_i = convert(Zonotope, Xk_i)
        sol_with_time[i] = ReachSet(Xk_i, t0, t1)
    end
    options = copy(sol.options)
    options[:n] += 1 # update state space dimension
    return ReachSolution(sol_with_time, options)
end

function project(sol::ReachSolution{Zonotope})
    N = length(sol.Xk)  # number of reach sets
    n = dim(first(sol.Xk).X) # state space dimension
    πsol = ReachSolution{Zonotope}(undef, N) # preallocated projected reachsets
    πvars = sol.options[:plot_vars] # variables for plotting
    @assert length(πvars) == 2

    if 0 ∈ πvars
        # add the time variable to the model (it is assumed it is not already
        # a variable in the model)
        sol = add_time(sol)
        n += 1
        πvars[first(indexin(0, πvars))] = n# time index is added in the end
    end

    M = sparse([1, 2], πvars, [1.0, 1.0], 2, n)
    for i in eachindex(sol.Xk)
        t0, t1 = sol.Xk[i].t0, sol.Xk[i].t1
        πsol_i = linear_map(M, sol.Xk[i].X)
        πsol_i = Zonotope(πsol_i.center, delete_zero_cols(πsol_i.generators))
        πsol_i = reduce_order(πsol_i, sol.options[:max_order])
        πsol[i] = ReachSet{Zonotope{Float64}, Float64}(πsol_i, t0, t1)
    end
    return πsol
end
