export solve

"""
    solve(problem::InitialValueProblem{ST, XT},
          opC::AbstractContinuousPost=default_continuous_post(ST),
          opD::Union{AbstractDiscretePost, Nothing}=nothing; kwargs...) where {ST, XT}

Solves a reachability problem for the  the given options.
If some options are not defined, we may fall back to default values.

### Input

- `system`    -- a (discrete or continuoues) system specification
- `options`   -- algorithm options for solving the problem
- `algorithm` -- (optional, default: dispatched on the system's type) the
                 reachability algorithm for the computation

### Output

A solution object whose content depends on the input options.

### Notes

To see all available input options, see
`keys(Reachability.available_keywords.dict)`.
"""
function solve(problem::InitialValueProblem{ST, XT},
               opC::AbstractContinuousPost=default_continuous_post(ST),
               opD::Union{AbstractDiscretePost, Nothing}=nothing; kwargs...) where {ST, XT}

    if iscontinuoussystem(ST)
        _solve_continuous(problem, opC; kwargs...)
    elseif ishybridsystem(ST)
        if opD == nothing
            opD = default_discrete_post(ST)
        end
        _solve_hybrid(problem, opC, opD, kwargs...)
    else
        error("a system of type $ST cannot be handled")
    end
end

function _solve_continuous(problem, op; kwargs...)
    @assert statedim(problem) == dim(problem.x0) "the state-space dimension should match the " *
    "dimension of the initial states, but they are of size $(statedim(p)) and $(problem.x0) respectively"

    # normalize system to canonical form if needed
    problem = IVP(normalize(problem.s), problem.x0)

    # initialize the algorithm-specific options
    #options = init!(op, problem, kwargs...)
    #options = Dict(kwargs...)
    T = kwargs[:T] # TODO check that key exists

    # run the continuous-post operator
    sol = post(op, problem, T)

    return sol
end

#=
function _solve_hybrid(problem, post, kwargs)
 # TODO
end
=#
