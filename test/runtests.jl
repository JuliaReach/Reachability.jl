#!/usr/bin/env julia
using LazySets, Reachability, MathematicalSystems

# fix namespace conflicts with MathematicalSystems
import LazySets.LinearMap

# compatibility between Julia versions
include("../src/compat.jl")
using Compat.Test

# in v0.7 and higher, the default Polyhedra library is used; in v0.6 CDDLib
# is used
@static if VERSION < v"0.7-"
    Pkg.add("CDDLib")
end

include("Systems/alltests.jl")
include("ReachSets/alltests.jl")
include("Reachability/alltests.jl")
include("Properties/alltests.jl")
