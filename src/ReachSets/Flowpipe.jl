"""
    Flowpipe{SN, RSN<:AbstractReachSet{SN}} <: AbstractSolution

Wrapper of a sequence of sets (the solution of a continuous-post operator).

### Fields

- `reachsets` -- the list of [`AbstractReachSet`](@ref)s
"""
struct Flowpipe{SN, RSN<:AbstractReachSet{SN}}
    reachsets::Vector{RSN}
end

# projection
function project(fp::Flowpipe, M::AbstractMatrix)
    return Flowpipe([project(X, M) for X in fp.reachsets])
end
