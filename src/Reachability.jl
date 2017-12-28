__precompile__()
"""
This is the main module and provides interfaces for specifying and solving reachability problems.
"""
module Reachability

include("logging.jl")
include("Systems/Systems.jl")
include("Utils/Utils.jl")
include("ReachSets/ReachSets.jl")
include("Properties/Properties.jl")
include("Transformations/Transformations.jl")

using Reexport, RecipesBase, Memento

@reexport using LazySets,
                Reachability.Utils,
                Reachability.Systems

using Reachability.ReachSets,
      Reachability.Properties,
      Reachability.Transformations

export LinearConstraintProperty, Clause,
       IntersectionProperty,
       Transformations

include("options.jl")
include("solve.jl")
include("plot_recipes.jl")

end # module
