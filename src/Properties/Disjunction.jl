"""
    Disjunction <: Property

Type that represents a disjunction of properties.

### Fields

- `disjuncts` -- vector of properties (elements are reordered by this type)
- `reorder`   -- flag to indicate whether shuffling is allowed

### Notes

The following formula characterizes whether a set ``X`` satisfies a disjunction
``𝑃 = 𝑃_1 ∨ 𝑃_2 ∨ … ∨ 𝑃_m``:

```math
    X \\models 𝑃 \\iff X \\models 𝑃_j \\text{ for some } 1 ≤ j ≤ m
```

If the `reorder` flag is set, the disjuncts may be reordered after each call to
[`check(𝑃::Disjunction, X::LazySet)`](@ref) as a heuristics to make subsequent
checks faster.
"""
struct Disjunction <: Property
    disjuncts::Vector{Property}
    reorder::Bool
end

# default constructor with activated reordering
Disjunction(disjuncts::Vector{<:Property}) = Disjunction(disjuncts, true)

"""
    check(𝑃::Disjunction, X::LazySet)::Bool

Check whether a convex set satisfies a disjunction of properties.

### Input

- `𝑃` -- disjunction of properties
- `X` -- convex set

### Output

`true` iff `X` satisfies the disjunction of properties `𝑃`.

### Notes

If the `𝑃.reorder` flag is set, the disjuncts may be reordered as a heuristics
to make subsequent checks faster.
Since we check satisfaction from left to right, we move the disjunct for which
satisfaction was established to the front.
"""
function check(𝑃::Disjunction, X::LazySet)::Bool
    for (i, conjunct) in enumerate(𝑃.disjuncts)
        if check(conjunct, X)
            _reorder!(𝑃, i)
            return true
        end
    end
    return false
end

function _reorder!(𝑃::Disjunction, i::Int)
    if !𝑃.reorder || i == 1
        return nothing
    end
    first = 𝑃.disjuncts[i]
    while i > 1
        𝑃.disjuncts[i] = 𝑃.disjuncts[i-1]
        i -= 1
    end
    𝑃.disjuncts[1] = first
    return nothing
end
