"""
    ReachSolution <: AbstractSolution

Solution of a reachability problem; more technically a wrapper of flowpipes
together with options.

### Fields

- `flowpipes` -- the list of [`Flowpipe`](@ref)s
- `options`   -- the dictionary of options
"""
struct ReachSolution <: AbstractSolution
    flowpipes::Vector{Flowpipe}
    options::AbstractOptions
end

# constructor with no options
ReachSolution(flowpipes::Vector{Flowpipe}) = ReachSolution(flowpipes, Options())

# constructor from a single flowpipe
ReachSolution(flowpipe::Flowpipe, options::Options) =
    ReachSolution([flowpipe], options)

# constructor from a single flowpipe no options
ReachSolution(flowpipe::Flowpipe) = ReachSolution([flowpipe], Options())

# projection
function project(rs::ReachSolution, M::AbstractMatrix)
    proj_flowpipes = [project(fp, M) for fp in rs.flowpipes]
    return ReachSolution(proj_flowpipes, rs.options)
end

# number of sets
function n_sets(rs::ReachSolution)
    return sum(fp -> length(fp.reachsets), rs.flowpipes)
end
