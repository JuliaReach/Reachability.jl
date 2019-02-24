"""
    SubsetProperty{N<:Real} <: Property

Type that represents a safety property characterized with a set of safe states.
A safety property is satisfied by a given set of states if the set of states is
fully contained in the set of safe states.

### Fields

- ``safe`` -- convex set representing the safe states
- ``witness`` -- witness point (empty vector if not set)
"""
mutable struct SubsetProperty{N<:Real} <: Property
    safe::LazySet
    witness::Vector{N}

    SubsetProperty{N}(safe::LazySet) where {N<:Real} = new(safe, N[])
end

# type-less convenience constructor
SubsetProperty(safe::LazySet{N}) where {N} = SubsetProperty{N}(safe)

"""
    check_property(set::LazySet, prop::SubsetProperty)::Bool

Checks whether a convex set satisfies a safety property.

### Input

- ``set``  -- convex set
- ``prop`` -- safety property with safe states

### Output

`true` iff the safety property is satisfied for the given set of states.
This is the case iff the set of states is a subset of the safe states.
"""
@inline function check_property(set::LazySet, prop::SubsetProperty)::Bool
    is_subset, witness = ⊆(set, prop.safe, true)
    if !is_subset
        # store violation witness
        prop.witness = witness
    end
    return is_subset
end
