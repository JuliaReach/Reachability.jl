"""
    ReachSolution{SN, RSN<:ReachSet{SN}} <: AbstractSolution

Type that wraps the solution of a reachability problem as a sequence of lazy
sets, and a dictionary of options.

### Fields

- `Xk`       -- the list of [`ReachSet`](@ref)s
- `options`  -- the dictionary of options
"""
struct ReachSolution{SN, RSN<:ReachSet{SN}} <: AbstractSolution
    Xk::Vector{RSN}
    options::AbstractOptions
end

# constructor with no options
ReachSolution(Xk::Vector{RSN}) where {SN, RSN<:ReachSet{SN}} =
    ReachSolution{SN, RSN<:ReachSet{SN}}(Xk, Options())

function project(rs::ReachSolution, M::AbstractMatrix)
    Yk = [project(X, M) for X in rs.Xk]
    return ReachSolution(Yk, rs.options)
end
