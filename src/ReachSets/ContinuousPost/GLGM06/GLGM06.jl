export GLGM06

struct GLGM06 <: ContinuousPost
    options::TwoLayerOptions

    function GLGM06(𝑂::Options)
        𝑂new = validate_and_wrap_options(𝑂, options_GLGM06())
        return new(𝑂new)
    end
end

# convenience constructor from pairs of symbols
GLGM06(𝑂::Pair{Symbol,<:Any}...) = GLGM06(Options(Dict{Symbol,Any}(𝑂)))

# default options (they are added in the function validate_and_wrap_options)
GLGM06() = GLGM06(Options())

include("init.jl")
include("post.jl")
include("reach.jl")
