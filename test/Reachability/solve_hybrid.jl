using Polyhedra # needed for some algorithms
using Reachability: constrained_dimensions

include("models/bouncing_ball.jl")
include("models/thermostat.jl")

models = [bouncing_ball, thermostat]

for model in models
    system, options = model()

    # --- reachability algorithms ---

    # default algorithm
    sol = solve(system, options)

    # concrete discrete-post operator
    sol = solve(system, options, BFFPSV18(), ConcreteDiscretePost())

    # lazy discrete-post operator
    sol = solve(system, options, BFFPSV18(), LazyDiscretePost())

    # overapproximating discrete-post operator
    sol = solve(system, options, BFFPSV18(), ApproximatingDiscretePost())

    # --- model analysis ---

    # constrained_dimensions
    HS = system.s
    if model == bouncing_ball
        @test constrained_dimensions(HS) == Dict{Int,Vector{Int}}(1=>[1, 2])
    elseif model == thermostat
        @test constrained_dimensions(HS) == Dict{Int,Vector{Int}}(1=>[1], 2=>[1])
    end
end
