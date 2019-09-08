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

function project(rs::ReachSet, M::AbstractMatrix)
    @assert dim(rs.X) == size(M, 2) "a projection of size $(size(M, 2)) is " *
        "incompatible with a set of dimension $(dim(rs.X))"
    return ReachSet(M * rs.X, rs.t_start, rs.t_end)
end

# =========================
# low-dimensional reach set
# =========================

struct SparseReachSet{SN}
    rs::ReachSet{SN}
    dimensions::Vector{Int}
end

function project(srs::SparseReachSet, M::AbstractMatrix)
    rs = srs.rs
    M_proj = M[:, srs.dimensions]
    return project(rs, M_proj)
end
