export TMJets

using IntervalArithmetic: IntervalBox

using StaticArrays: SVector

struct TMJets <: ContinuousPost
    options::TwoLayerOptions

    function TMJets(𝑂::Options)
        𝑂new = validate_and_wrap_options(𝑂, options_TMJets())
        return new(𝑂new)
    end
end

# convenience constructor from pairs of symbols
TMJets(𝑂::Pair{Symbol, <:Any}...) = TMJets(Options(Dict{Symbol, Any}(𝑂)))

# default options (they are added in the function validate_and_wrap_options)
TMJets() = TMJets(Options())

include("init.jl")
include("post.jl")
include("project.jl")
include("reach.jl")
