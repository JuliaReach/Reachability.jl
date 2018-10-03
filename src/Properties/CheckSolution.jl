"""
    CheckSolution

Type that wraps the solution of a property checking problem.

### Fields

- `satisfied` -- is the property satisfied?
- `violation` -- step at which the property is violated (-1 otherwise)
- `options`   -- the dictionary of options

### Notes

This type contains the answer if the property is satisfied, and if not, it
contains the index at which the property might be violated for the first time.
"""
struct CheckSolution <: AbstractSolution
    satisfied::Bool
    violation::Int
    options::Options
end

# constructor with no options
CheckSolution(satisfied::Bool, violation::Int) =
    CheckSolution(satisfied, violation, Options())
