"""
Abstract supertype of properties that can be checked.

Every concrete subtype should provide the following function:
  - `check(𝑃::Property, X::LazySet)::Bool`
"""
abstract type Property end
