using Test, LazySets, Reachability, MathematicalSystems, MathematicalPredicates,
      LinearAlgebra, SparseArrays, IntervalMatrices

# fix namespace conflicts with MathematicalSystems
using LazySets: LinearMap, AffineMap, ResetMap, Interval

include("Systems/alltests.jl")
include("ReachSets/alltests.jl")
include("Reachability/alltests.jl")
include("Properties/alltests.jl")
