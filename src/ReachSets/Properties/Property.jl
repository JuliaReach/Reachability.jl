"""
Abstract supertype of properties that can be checked.

Every concrete subtype should provide the following functions:
  - `inout_map_property(::Property, ::AbstractVector{<:AbstractVector{Int}},
                        ::AbstractVector{Int}, ::Int)::Property`
  - `check_property(::LazySet, ::Property)::Bool`
"""
abstract type Property end
