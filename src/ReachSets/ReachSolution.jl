"""
    ReachSolution{S<:LazySet} <: AbstractSolution

Type that wraps the solution of a reachability problem as a sequence of lazy
sets, and a dictionary of options.

### Fields

- `Xk`       -- the list of [`ReachSet`](@ref)s
- `options`  -- the dictionary of options
"""
struct ReachSolution{S<:LazySet} <: AbstractSolution
    Xk::Vector{<:ReachSet{S}}
    options::Options
end

# constructor with no options
ReachSolution(Xk::Vector{<:ReachSet{S}}) where {S<:LazySet} =
    ReachSolution{S}(Xk, Options())
