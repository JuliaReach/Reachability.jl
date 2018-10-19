#!/usr/bin/env julia
using LazySets, Reachability, MathematicalSystems

# compatibility between Julia versions
include("../src/compat.jl")
using Compat.Test

include("Systems/alltests.jl")
include("ReachSets/alltests.jl")
include("Reachability/alltests.jl")
include("Properties/alltests.jl")
