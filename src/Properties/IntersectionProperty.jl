"""
    IntersectionProperty{N<:Real} <: Property

Type that represents a safety property characterized with a set of bad states.
A safety property is satisfied by a given set of states if the intersection with
the set of bad states is empty.

### Fields

- ``bad`` -- convex set representing the bad states
- ``witness`` -- witness point (empty vector if not set)
"""
mutable struct IntersectionProperty{N<:Real} <: Property
    bad::LazySet
    witness::Vector{N}

    IntersectionProperty{N}(bad::LazySet) where {N<:Real} = new(bad, N[])
end

# type-less convenience constructor
IntersectionProperty(bad::LazySet{N}) where {N} = IntersectionProperty{N}(bad)

"""
    inout_map_property(prop::IntersectionProperty,
                       partition::AbstractVector{<:AbstractVector{Int}},
                       blocks::AbstractVector{Int},
                       n::Int
                      )::IntersectionProperty

Map an `IntersectionProperty` to the dimensions of analyzed blocks.

### Input

- `prop`      -- property
- `partition` -- block partition; elements are start and end indices of a block
- `blocks`    -- list of all output block indices in the partition
- `n`         -- total number of input dimensions

### Output

A new property of reduced dimension.

### Notes

If the dimension is not reduced, we keep the original set.
Otherwise, the dimension reduction is achieved with a `LinearMap`.
"""
function inout_map_property(prop::IntersectionProperty{N},
                            partition::AbstractVector{<:AbstractVector{Int}},
                            blocks::AbstractVector{Int},
                            n::Int
                           )::IntersectionProperty{N} where {N<:Real}
    proj = projection_map(partition, blocks)
    if length(proj) == n
        # no change in the dimension, copy the old property (keep the set)
        return IntersectionProperty(prop.bad)
    else
        M = sparse(proj, proj, ones(N, length(proj)), n, n)
        return IntersectionProperty(M * prop.bad)
    end
end

"""
    check_property(set::LazySet, prop::IntersectionProperty)::Bool

Checks whether a convex set satisfies a safety property.

### Input

- ``set``  -- convex set
- ``prop`` -- safety property with bad states

### Output

`true` iff the safety property is satisfied for the given set of states.
This is the case iff the set of states does not intersect with the bad states.
"""
@inline function check_property(set::LazySet, prop::IntersectionProperty)::Bool
    nonempty_intersection, witness = is_intersection_empty(set, prop.bad, true)
    if nonempty_intersection
        # store violation witness
        prop.witness = witness
    end
    return !nonempty_intersection
end
