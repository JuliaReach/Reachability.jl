__precompile__()
"""
Module for defining and checking properties.
"""
module Properties

using LazySets

# ==============================
# Property struct and evaluation
# ==============================
include("Property.jl")
export Property,
       check

include("Conjunction.jl")
include("Disjunction.jl")
export Conjunction,
       Disjunction

include("IntersectionProperty.jl")
export IntersectionProperty

include("SubsetProperty.jl")
export SubsetProperty

end # module
