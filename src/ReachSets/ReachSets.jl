__precompile__()
"""
Module with reachability algorithms.
"""
module ReachSets

using ..Utils
using LazySets, MathematicalSystems, MathematicalPredicates, HybridSystems,
      Expokit, Optim, LinearAlgebra, SparseArrays
import IntervalMatrices
using IntervalMatrices: IntervalMatrix

import ProgressMeter
using ProgressMeter: Progress
# no-op progress meter for unbounded time
function ProgressMeter.update!(::Nothing, ::Int) end

# fix namespace conflicts
using LazySets: LinearMap, AffineMap, ResetMap, Interval
using LazySets.Approximations: AbstractDirections

import LazySets: set
import LazySets.Approximations: project

using Reachability: info, warn,
                    @timing,
                    Options, OptionSpec, TwoLayerOptions, AbstractOptions,
                    validate_and_wrap_options, print_option_spec,
                    haskey_specified,
                    validate_solver_options_and_add_default_values!

# ========================================
# Discretize and compute bloating factors
# ========================================
include("discretize.jl")

# ==========================
# Property checking results
# ==========================
include("CheckSolution.jl")

export CheckSolution

# =============================
# Property checking algorithms
# =============================

# dictionary of registered algorithms
available_algorithms_check = Dict{String, Dict{String, Any}}()

include("ContinuousPost/BFFPSV18/check_blocks.jl")
include("ContinuousPost/BFFPSV18/check_property.jl")
include("ContinuousPost/BFFPSV18/partitions.jl")
include("ContinuousPost/BFFPSV18/compute_dimensions.jl")

# "explicit" backends
push!(available_algorithms_check, "explicit_blocks"=>Dict("func"=>check_blocks,
                                                          "is_explicit"=>true))

export available_algorithms_check,
       check_property

# ====================================================
# Algorithms to find a threshold for property checking
# ====================================================
include("tune.jl")
export tune_δ

# =====================
# Reachability results
# =====================
include("ReachSet.jl")
include("Flowpipe.jl")
include("ReachSolution.jl")

export AbstractReachSet, ReachSet, SparseReachSet,
       set, time_start, time_end,
       Flowpipe,
       ReachSolution,
       n_sets

# ===============
# Post operators
# ===============
include("PostOperator.jl")
include("ContinuousPost/ContinuousPost.jl")
include("DiscretePost/DiscretePost.jl")

export PostOperator,
       ContinuousPost,
       DiscretePost,
       init, init!,
       post,
       tube⋂inv

# ==========================
# Continuous post operators
# ==========================
include("ContinuousPost/BFFPSV18/BFFPSV18.jl")
include("ContinuousPost/BFFPSV18/reach.jl")
include("ContinuousPost/BFFPSV18/reach_blocks.jl")
include("ContinuousPost/BFFPSV18/reach_blocks_wrapping_effect.jl")

include("ContinuousPost/GLGM06/GLGM06.jl")

include("ContinuousPost/TMJets/TMJets.jl")

include("ContinuousPost/BFFPS19/BFFPS19.jl")

include("ContinuousPost/ASB07/ASB07.jl")

include("ContinuousPost/ASB07_decomposed/ASB07_decomposed.jl")

# ========================
# Reachability Algorithms
# ========================
import Reachability.check_aliases_and_add_default_value!

# dictionary of registered algorithms
available_algorithms = Dict{String, Dict{String, Any}}()

# "explicit" backends
push!(available_algorithms, "explicit_blocks"=>Dict("func"=>reach_blocks!,
                                                    "is_explicit"=>true))

push!(available_algorithms,
      "wrap"=>Dict("func"=>reach_blocks_wrapping_effect!,
                   "is_explicit"=>true))

export available_algorithms

# ==========================
# Discrete post operators
# ==========================
include("DiscretePost/LazyDiscretePost.jl")
include("DiscretePost/ConcreteDiscretePost.jl")
include("DiscretePost/DecomposedDiscretePost.jl")

# ==============================================
# Projection of the reach set in two dimensions
# ==============================================
include("project_reach.jl")

export project_reach

end # module ReachSets
