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

include("LinearConstraintProperty.jl")
export LinearConstraintProperty,
       Clause

include("IntersectionProperty.jl")
export IntersectionProperty

include("SubsetProperty.jl")
export SubsetProperty

end # module
