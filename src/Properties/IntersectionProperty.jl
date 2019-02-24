"""
    IntersectionProperty{N<:Real} <: Property

Type that represents a safety property characterized by a set of bad states.
The property is satisfied by a given set of states if the intersection with the
set of bad states is empty.

### Fields

- `bad`     -- convex set representing the bad states
- `witness` -- witness point (empty vector if not set)
"""
mutable struct IntersectionProperty{N<:Real} <: Property
    bad::LazySet
    witness::Vector{N}

    IntersectionProperty{N}(bad::LazySet) where {N<:Real} = new(bad, N[])
end

# type-less convenience constructor
IntersectionProperty(bad::LazySet{N}) where {N<:Real} =
    IntersectionProperty{N}(bad)

"""
    check(𝑃::IntersectionProperty, X::LazySet)::Bool

Checks whether a convex set is disjoint from the set of bad states.

### Input

- `𝑃` -- safety property with bad states
- `X` -- convex set

### Output

`true` iff the given set of states does not intersect with the bad states.
"""
@inline function check(𝑃::IntersectionProperty, X::LazySet)::Bool
    empty_intersection, witness = is_intersection_empty(X, 𝑃.bad, true)
    if !empty_intersection
        # store violation witness
        𝑃.witness = witness
    end
    return empty_intersection
end
