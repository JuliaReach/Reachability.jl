abstract type AbstractReachSet{SN} end

"""
    ReachSet{SN}

Type that wraps a reach set, which is a set representing (an approximation of)
the reachable states for a given time interval.

### Fields

- `X`       -- set
- `t_start` -- time interval lower bound
- `t_end`   -- time interval upper bound
"""
struct ReachSet{SN} <: AbstractReachSet{SN}
    X::SN
    t_start::Float64
    t_end::Float64
end

function set(rs::ReachSet)
    return rs.X
end

function time_start(rs::ReachSet)
    return rs.t_start
end

function time_end(rs::ReachSet)
    return rs.t_end
end

function substitute(rs::ReachSet; set=set(rs), time_start=time_start(rs),
                    time_end=time_end(rs))
    return ReachSet(set, time_start, time_end)
end

function project(rs::ReachSet, M::AbstractMatrix)
    @assert dim(rs.X) == size(M, 2) "a projection of size $(size(M, 2)) is " *
        "incompatible with a set of dimension $(dim(rs.X))"
    return ReachSet(M * rs.X, rs.t_start, rs.t_end)
end

# =========================
# low-dimensional reach set
# =========================

struct SparseReachSet{SN} <: AbstractReachSet{SN}
    rs::ReachSet{SN}
    dimensions::AbstractVector{Int}
end

# flattened constructor
function SparseReachSet(X, t_start::Float64, t_end::Float64,
                        dimensions::AbstractVector{Int})
    return SparseReachSet(ReachSet(X, t_start, t_end), dimensions)
end

function set(srs::SparseReachSet)
    return set(srs.rs)
end

function time_start(srs::SparseReachSet)
    return time_start(srs.rs)
end

function time_end(srs::SparseReachSet)
    return time_end(srs.rs)
end

function substitute(srs::SparseReachSet; set=set(srs),
                    time_start=time_start(srs), time_end=time_end(srs))
    substituted_rs = substitute(srs.rs, set=set, time_start=time_start, time_end=time_end)
    return SparseReachSet(substituted_rs, srs.dimensions)
end

function project(srs::SparseReachSet, M::AbstractMatrix)
    rs = srs.rs
    M_proj = M[:, srs.dimensions]
    return project(rs, M_proj)
end
