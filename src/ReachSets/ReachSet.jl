"""
    ReachSet{SN}

Type that wraps a reach set, which is a set representing (an approximation of)
the reachable states for a given time interval.

### Fields

- `X`       -- set
- `t_start` -- time interval lower bound
- `t_end`   -- time interval upper bound
"""
struct ReachSet{SN}
    X::SN
    t_start::Float64
    t_end::Float64
end
