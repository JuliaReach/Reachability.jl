export ASB07_decomposed

struct ASB07_decomposed <: AbstractContinuousPost
    options::TwoLayerOptions

    function ASB07_decomposed(𝑂::Options)
        𝑂new = validate_and_wrap_options(𝑂, options_ASB07_decomposed())
        return new(𝑂new)
    end
end

# convenience constructor from pairs of symbols
ASB07_decomposed(𝑂::Pair{Symbol,<:Any}...) =
    ASB07_decomposed(Options(Dict{Symbol,Any}(𝑂)))

# default options (they are added in the function validate_and_wrap_options)
ASB07_decomposed() = ASB07_decomposed(Options())

include("options.jl")
include("init.jl")
include("post.jl")
include("reach.jl")
