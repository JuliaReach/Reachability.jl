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
    check(𝑃::Conjunction, X::LazySet{N}; witness::Bool=false) where {N<:Real}

Check whether a convex set satisfies a conjunction of properties.

### Input

- `𝑃`       -- conjunction of properties
- `X`       -- convex set
- `witness` -- (optional, default: `false`) flag for returning a counterexample
               if the property is violated

### Output

* If `witness` option is deactivated: `true` iff `X` satisfies the property `𝑃`
* If `witness` option is activated:
  * `(true, [])` iff `X` satisfies the property `𝑃`
  * `(false, v)` iff `X` does not satisfy the property `𝑃` with witness `v`

### Notes

By convention, the empty conjunction is equivalent to `true` and hence is
satisfied by any set.
"""
function check(𝑃::Conjunction, X::LazySet{N};
               witness::Bool=false) where {N<:Real}
    for conjunct in 𝑃.conjuncts
        result = check(conjunct, X; witness=witness)
        if (witness && !result[1]) || !result
            return result
        end
    end
    return witness ? (true, N[]) : true
end
