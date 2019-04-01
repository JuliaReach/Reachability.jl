export ASB08

# load packages
using LazySets.Approximations
using Reachability.ReachSets: Φ₁
using IntervalArithmetic
using LazySets: Interval, translate
using TaylorSeries
using TaylorSeries: gradient, jacobian, hessian, derivative
const ∂ = derivative

struct ASB08 <: ContinuousPost
    options::TwoLayerOptions

    function ASB08(𝑂::Options)
        𝑂new = validate_and_wrap_options(𝑂, options_ASB08())
        return new(𝑂new)
    end
end

# convenience constructor from pairs of symbols
ASB08(𝑂::Pair{Symbol,<:Any}...) = ASB08(Options(Dict{Symbol,Any}(𝑂)))

# default options (they are added in the function validate_and_wrap_options)
ASB08() = ASB08(Options())

include("init.jl")
include("post.jl")
include("reach.jl")
