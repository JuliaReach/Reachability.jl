__precompile__()
"""
Reachability analysis of high-dimensional affine systems with nondeterministic inputs.

## Introduction

This module provides an interface to set-based reachability algorithms for affine
systems with nondeterministic inputs. It implements several reachability
algorithms using decomposition into low-dimensional subsystems. The code is optimized
for high-dimensional linear systems.

ReachSets.jl relies in the following modules:

|Module name|Features|
|-----------|--------|
|`LazySets.jl`|Set-based representations and operations in the lazy framework. Support vector and support function calculus.|
|`MathematicalSystems.jl`|Representation and transformation of affine systems with non-deterministic inputs.|

The following dependencies are optional and provide additional functionality:

|Optional module name | Features|
|------|------|
|`Expokit.jl` | Fast computation of matrix exponential of very large and very sparse matrices.|
|`Polyhedra.jl`, `CDDLib.jl` | Julia Polyhedra library and the CDD backend.|

The visualization routines are deferred to the corresponding `LazySets.jl`, such as
`plot_polygon`, implemented in `Polygon.jl`.

This module is accompanied with a set of unit tests, available in the `test` folder.

## Available algorithms

The dictionary ``available_algorithms`` specifies the registered algorithms. Each
algorithm is registered with its name, and a dictionary value which any number
of entries. Most common entries are:

|Key|Description|
|-----|----|
|``'func'``|the corresponding function's handle|
|``'is_explicit'`` |if true, the output contains explicitly computed polygons; otherwise, the lazy-set framework is returned|

The algorithms implemented in `ReachSets/` are summarized in the following
table.

|Function name|Description|
|-------------|------|
|``lazy``     | Lazy computation (all blocks)|
|``explicit`` | Explicit computation|

For additional information on other algorithms see the documentation in their
corresponding source files.

## Projections

The final goal of any reachability computation is a sequence of points (vertices)
describing the evolution of one or more variables in time, as a function of
other variables, or a linear combination of them. This information is required
for plotting the output in two or three dimensions.

This module provides the function ``project_reach`` for the purpose of projecting
a high-dimensional explicit or lazy set computation into a low number of variables
efficiently. Additionally, the time variable can be passed (use `0`), and an arbitrary
linear combination of state variables can be handled as well.

## Helper functions and macros

This module provides commonly used

|Macro or function name| Description|
|-----|-----|
|``@filename_to_png``|Convert the current script name's suffix to ".png"|
|``@block_id v``| Return the block number associated to a given variable.|


AUTHORS:

- Marcelo Forets
- Christian Schilling
- Frederic Viry
"""
module ReachSets

using ..Utils
using LazySets, MathematicalSystems, HybridSystems
using Expokit, Optim, ProgressMeter

# fix namespace conflicts with MathematicalSystems
using LazySets: LinearMap
using Reachability: info

include("../compat.jl")

using LazySets.Approximations: symmetric_interval_hull,
                               decompose,
                               overapproximate,
                               box_approximation
using Reachability: @timing,
                    Options, OptionSpec, TwoLayerOptions,
                    validate_and_wrap_options, print_option_spec,
                    validate_solver_options_and_add_default_values!

# ========================================
# Discretize and compute bloating factors
# ========================================
include("discretize.jl")

export discretize

# ==============================
# Property struct and evaluation
# ==============================
include("Properties/Property.jl")
export Property,
       inout_map_property

include("Properties/LinearConstraintProperty.jl")
export LinearConstraintProperty,
       Clause

include("Properties/IntersectionProperty.jl")
export IntersectionProperty

include("Properties/SubsetProperty.jl")
export SubsetProperty

# ==========================
# Property checking results
# ==========================
include("Properties/CheckSolution.jl")

export CheckSolution

# =============================
# Property checking algorithms
# =============================

# dictionary of registered algorithms
available_algorithms_check = Dict{String, Dict{String, Any}}()

include("ContinuousPost/BFFPSV18/check_blocks.jl")
include("ContinuousPost/BFFPSV18/check_property.jl")
include("ContinuousPost/BFFPSV18/partitions.jl")
include("ContinuousPost/BFFPSV18/inout_map_reach.jl")

# "explicit" backends
push!(available_algorithms_check, "explicit_blocks"=>Dict("func"=>check_blocks,
                                                          "is_explicit"=>true))

export available_algorithms_check,
       check_property

# ====================================================
# Algorithms to find a threshold for property checking
# ====================================================
include("Properties/tune.jl")
export tune_δ

# =====================
# Reachability results
# =====================
include("ReachSet.jl")
include("ReachSolution.jl")

export ReachSet,
       ReachSolution

# ===============
# Post operators
# ===============
include("PostOperator.jl")
include("ContinuousPost/ContinuousPost.jl")
include("DiscretePost/DiscretePost.jl")

export PostOperator,
       ContinuousPost,
       DiscretePost,
       init,
       post,
       tube⋂inv!

# ==========================
# Continuous post operators
# ==========================
include("ContinuousPost/BFFPSV18/BFFPSV18.jl")
include("ContinuousPost/BFFPSV18/reach.jl")
include("ContinuousPost/BFFPSV18/reach_blocks.jl")
include("ContinuousPost/BFFPSV18/reach_blocks_wrapping_effect.jl")

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

# ==============================================
# Projection of the reach set in two dimensions
# ==============================================
include("project_reach.jl")

export project_reach,
       project

end # module ReachSets
