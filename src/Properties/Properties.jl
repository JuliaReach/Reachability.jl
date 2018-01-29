__precompile__()
"""
This module checks a given safety property without explicitly computing the reachable states.
For example, to check whether x_42 is less than 5, we just need to track the maximum values of x_42, i.e., evaluate the support of x_42 in the positive infinity direction.
This is cheaper than computing the reachable states also for the negative direction (and also for the other variable in the respective block).
"""
module Properties

import LazySets.Approximations:decompose,
                               overapproximate

using LazySets, ..Systems, Expokit

# ==============================
# Property struct and evaluation
# ==============================
include("Property.jl")
include("LinearConstraintProperty.jl")
include("IntersectionProperty.jl")
include("SubsetProperty.jl")

# dictionary of registered algorithms
available_algorithms = Dict{String, Dict{String, Any}}()

# "explicit" backends
include("check_explicit_block.jl")
push!(available_algorithms, "explicit_block"=>Dict("func"=>check_explicit_block!,
                                                   "is_full"=>false,
                                                   "is_explicit"=>true))

include("check_explicit_blocks.jl")
push!(available_algorithms, "explicit_blocks"=>Dict("func"=>check_explicit_blocks!,
                                                    "is_full"=>false,
                                                    "is_explicit"=>true))

include("check_property.jl")

export check_property

end #module Properties