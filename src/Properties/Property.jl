export Property

"""
Abstract supertype of properties that can be checked.

Every concrete subtype should provide the following function:
  - `check_property(::LazySet, ::Property)::Bool`
"""
abstract type Property end
