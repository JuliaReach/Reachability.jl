# out-of-place initialization
init(𝒜::GLGM06, 𝑆::AbstractSystem, 𝑂::Options) = init!(𝒜, 𝑆, copy(𝑂))

function options_GLGM06()

    𝑂spec = Vector{OptionSpec}()

    # step size
    push!(𝑂spec, OptionSpec(:δ, 1e-2, domain=Float64, aliases=[:sampling_time],
                            domain_check=(v  ->  v > 0.), info="time step"))

    # discretization
    push!(𝑂spec, OptionSpec(:discretization, "forward", domain=String,
                            info="model for bloating/continuous time analysis"))

    push!(𝑂spec, OptionSpec(:sih_method, "concrete", domain=String,
                            info="method to compute the symmetric interval hull in discretization"))

    push!(𝑂spec, OptionSpec(:exp_method, "base", domain=String,
                            info="method to compute the matrix exponential"))

    # approximation options
    push!(𝑂spec, OptionSpec(:max_order, 10, domain=Int,
                            info="maximum allowed order of zonotopes"))

    return 𝑂spec
end

# in-place initialization
function init!(::GLGM06, 𝑆::AbstractSystem, 𝑂::Options)

    # state dimension
    𝑂[:n] = statedim(𝑆)

    # adds default values for unspecified options
    𝑂init = validate_solver_options_and_add_default_values!(𝑂)

    return 𝑂init
end
