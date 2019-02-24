"""
    Conjunction <: Property

Type that represents a conjunction of properties.

### Fields

- `conjuncts` -- vector of properties

### Notes

The following formula characterizes whether a set ``X`` satisfies a disjunction
``𝑃 = 𝑃_1 ∧ 𝑃_2 ∧ … ∧ 𝑃_m``:

```math
    X \\models 𝑃 \\iff X \\models 𝑃_j \\text{ for all } 1 ≤ j ≤ m
```
"""
struct Conjunction <: Property
    conjuncts::Vector{Property}
end

"""
    check(𝑃::Conjunction, X::LazySet)::Bool

Check whether a convex set satisfies a conjunction of properties.

### Input

- `𝑃` -- conjunction of properties
- `X` -- convex set

### Output

`true` iff `X` satisfies the conjunction of properties `𝑃`.
"""
function check(𝑃::Conjunction, X::LazySet)::Bool
    for conjunct in 𝑃.conjuncts
        if !check(conjunct, X)
            return false
        end
    end
    return true
end
