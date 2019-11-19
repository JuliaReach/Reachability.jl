"""
    AbstractPostOperator

Abstract supertype of all post operators.

### Notes

All post operators should provide the following methods:
```julia
init(op::AbstractPostOperator, system, options)

post(op::AbstractPostOperator, system, ...)
```
"""
abstract type AbstractPostOperator end
