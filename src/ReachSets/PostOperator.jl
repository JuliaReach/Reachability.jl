"""
    PostOperator

Abstract supertype of all post operators.

### Notes

All post operators should provide the following methods:
```julia
init(op::PostOperator, system, options)

post(op::PostOperator, system, ...)
```
"""
abstract type PostOperator end
