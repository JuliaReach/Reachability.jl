"""
    ReachSet{S<:LazySet, N<:Real}

Type that wraps a reach set, which is a set representing (an approximation of)
the reachable states for a given time interval.

### Fields

- `X`       -- set
- `t_start` -- time interval lower bound
- `t_end`   -- time interval upper bound
- `k`   -- time step
"""
struct ReachSet{S<:LazySet, N<:Real}
    X::S
    t_start::N
    t_end::N
    k::Int
end
