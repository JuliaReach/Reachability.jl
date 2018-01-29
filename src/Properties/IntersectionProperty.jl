export IntersectionProperty

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
# convenience constructor for Float64
IntersectionProperty(bad::LazySet) = IntersectionProperty{Float64}(bad)

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
