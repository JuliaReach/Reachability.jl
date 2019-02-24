"""
    BadStatesProperty{N<:Real} <: Property

Type that represents a safety property characterized by a set of bad states.
The property is satisfied by a given set of states if the intersection with the
set of bad states is empty.

### Fields

- `bad`     -- convex set representing the bad states
- `witness` -- witness point (empty vector if not set)

### Notes

The following formula characterizes whether a set ``X`` satisfies a safety
property characterized by a set of bad states 𝑃:

```math
    X \\models 𝑃 \\iff X ∩ 𝑃.\\texttt{bad} = ∅
```
"""
mutable struct BadStatesProperty{N<:Real} <: Property
    bad::LazySet
    witness::Vector{N}

    BadStatesProperty{N}(bad::LazySet) where {N<:Real} = new(bad, N[])
end

# type-less convenience constructor
BadStatesProperty(bad::LazySet{N}) where {N<:Real} =
    BadStatesProperty{N}(bad)

"""
    check(𝑃::BadStatesProperty, X::LazySet)::Bool

Checks whether a convex set is disjoint from the set of bad states.

### Input

- `𝑃` -- safety property with bad states
- `X` -- convex set

### Output

`true` iff the given set of states does not intersect with the set of bad
states.
"""
@inline function check(𝑃::BadStatesProperty, X::LazySet)::Bool
    empty_intersection, witness = is_intersection_empty(X, 𝑃.bad, true)
    if !empty_intersection
        # store violation witness
        𝑃.witness = witness
    end
    return empty_intersection
end
