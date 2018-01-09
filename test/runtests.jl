#!/usr/bin/env julia
include("../src/Reachability.jl")
using Base.Test, LazySets, Reachability

include("Systems/alltests.jl")
include("ReachSets/alltests.jl")
include("Reachability/alltests.jl")
include("Properties/alltests.jl")
