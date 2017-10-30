__precompile__()
"""
This is the main module and provides interfaces for specifying and solving reachability problems.
"""
module Reachability

include("Systems/Systems.jl")
include("Utils/Utils.jl")
include("ReachSets/ReachSets.jl")
include("Properties/Properties.jl")
include("Transformations/Transformations.jl")

using Reexport, RecipesBase

@reexport using LazySets,
                Reachability.Utils,
                Reachability.Systems

using Reachability.ReachSets,
      Reachability.Properties,
      Reachability.Transformations

export Property, Clause

include("options.jl")
include("solve.jl")
include("plot_recipes.jl")

end # module
