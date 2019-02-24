#!/usr/bin/env julia
using LazySets, Reachability, MathematicalSystems,
      LinearAlgebra, SparseArrays

# fix namespace conflicts with MathematicalSystems
import LazySets.LinearMap

# compatibility between Julia versions
include("../src/compat.jl")
using Compat.Test

# Note: CDDLib is only required for v0.6; for v0.7 or later the default Polyhedra
# library is used and CDDLib can be removed from REQUIRE

include("Systems/alltests.jl")
include("ReachSets/alltests.jl")
include("Reachability/alltests.jl")
include("Properties/alltests.jl")
