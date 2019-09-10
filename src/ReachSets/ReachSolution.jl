"""
    ReachSolution{SN, RSN<:AbstractReachSet{SN}} <: AbstractSolution

Type that wraps the solution of a reachability problem as a sequence of lazy
sets, and a dictionary of options.

### Fields

- `Xk`       -- the list of [`AbstractReachSet`](@ref)s
- `options`  -- the dictionary of options
"""
struct ReachSolution{SN, RSN<:AbstractReachSet{SN}} <: AbstractSolution
    Xk::Vector{RSN}
    options::AbstractOptions
end

# constructor with no options
ReachSolution(Xk::Vector{RSN}) where {SN, RSN<:AbstractReachSet{SN}} =
    ReachSolution{SN, RSN}(Xk, Options())

function project(rs::ReachSolution, M::AbstractMatrix)
    Yk = [project(X, M) for X in rs.Xk]
    return ReachSolution(Yk, rs.options)
end
