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

include("BadStatesProperty.jl")
export BadStatesProperty

include("SafeStatesProperty.jl")
export SafeStatesProperty

end # module
