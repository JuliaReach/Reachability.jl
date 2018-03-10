__precompile__()
"""
This is the main module and provides interfaces for specifying and solving reachability problems.
"""
module Reachability

using Reexport, RecipesBase, Memento
@reexport using LazySets

include("logging.jl")
include("systems.jl")
include("Utils/Utils.jl")
include("ReachSets/ReachSets.jl")
include("Properties/Properties.jl")
include("Transformations/Transformations.jl")

@reexport using Reachability.Utils

using Reachability.ReachSets,
      Reachability.Properties,
      Reachability.Transformations

export Properties, LinearConstraintProperty, Clause,
       IntersectionProperty,
       SubsetProperty,
       Transformations

include("options.jl")
include("solve.jl")
include("plot_recipes.jl")

end # module
