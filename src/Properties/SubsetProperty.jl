"""
    SubsetProperty{N<:Real} <: Property

Type that represents a safety property characterized by a set of safe states.
The property is satisfied by a given set of states ``X`` if ``X`` is fully
contained in the set of safe states.

### Fields

- `safe`    -- convex set representing the safe states
- `witness` -- witness point (empty vector if not set)
"""
mutable struct SubsetProperty{N<:Real} <: Property
    safe::LazySet
    witness::Vector{N}

    SubsetProperty{N}(safe::LazySet) where {N<:Real} = new(safe, N[])
end

# type-less convenience constructor
SubsetProperty(safe::LazySet{N}) where {N<:Real} = SubsetProperty{N}(safe)

"""
    check(𝑃::SubsetProperty, X::LazySet)::Bool

Checks whether a convex set is contained in the set of safe states.

### Input

- `𝑃` -- safety property with safe states
- `X` -- convex set

### Output

`true` iff the given set of states is a subset of the set of safe states.
"""
@inline function check(𝑃::SubsetProperty, X::LazySet)::Bool
    is_subset, witness = ⊆(X, 𝑃.safe, true)
    if !is_subset
        # store violation witness
        𝑃.witness = witness
    end
    return is_subset
end
