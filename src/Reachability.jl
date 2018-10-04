__precompile__()
"""
This is the main module and provides interfaces for specifying and solving reachability problems.
"""
module Reachability

using Reexport, RecipesBase, Memento, MathematicalSystems, HybridSystems
@reexport using LazySets

include("logging.jl")
include("Utils/Utils.jl")
include("options.jl")
include("ReachSets/ReachSets.jl")
include("Transformations/Transformations.jl")

@reexport using Reachability.Utils

using Reachability.ReachSets,
      Reachability.Transformations

export project,
       LinearConstraintProperty, Clause,
       IntersectionProperty,
       SubsetProperty,
       Transformations

include("solve.jl")
include("plot_recipes.jl")

end # module
