__precompile__()
"""
This module checks a given safety property without explicitly computing the reachable states.
For example, to check whether x_42 is less than 5, we just need to track the maximum values of x_42, i.e., evaluate the support of x_42 in the positive infinity direction.
This is cheaper than computing the reachable states also for the negative direction (and also for the other variable in the respective block).
"""
module Properties

import LazySets.Approximations:decompose,
                               overapproximate
import Reachability.tocc

using LazySets, MathematicalSystems, ..Utils, Expokit, ProgressMeter

# ==============================
# Property struct and evaluation
# ==============================
include("Property.jl")
export Property,
       inout_map_property

include("LinearConstraintProperty.jl")
export LinearConstraintProperty,
       Clause

include("IntersectionProperty.jl")
export IntersectionProperty

include("SubsetProperty.jl")
export SubsetProperty

# dictionary of registered algorithms
available_algorithms = Dict{String, Dict{String, Any}}()

# "explicit" backends
include("check_blocks.jl")
push!(available_algorithms, "explicit_blocks"=>Dict("func"=>check_blocks,
                                                    "is_full"=>false,
                                                    "is_explicit"=>true))

include("check_property.jl")
export check_property

# ====================================================
# Algorithms to find a threshold for property checking
# ====================================================
include("tune.jl")
export tune_δ

end #module Properties
