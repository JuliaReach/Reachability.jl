export GLGM06

#=
struct GLGM06 <: AbstractContinuousPost
    options::TwoLayerOptions

    function GLGM06(𝑂::Options)
        𝑂new = validate_and_wrap_options(𝑂, options_GLGM06())
        return new(𝑂new)
    end
end
=#

@with_kw struct GLGM06 <: AbstractContinuousPost
    δ::Float64 = 1e-2
    discretization::String = "forward"
    sih_method::String = "concrete"
    exp_method::String = "base"
    max_order::Int = 10
end

# convenience constructor from pairs of symbols
#GLGM06(𝑂::Pair{Symbol,<:Any}...) = GLGM06(Options(Dict{Symbol,Any}(𝑂)))

# default options (they are added in the function validate_and_wrap_options)
#GLGM06() = GLGM06(Options())

include("init.jl")
include("post.jl")
include("reach.jl")
include("project.jl")
