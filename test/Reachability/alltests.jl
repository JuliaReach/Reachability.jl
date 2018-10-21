@time @testset "Reachability.options" begin include("unit_options.jl") end
@time @testset "Reachability.solve_continuous" begin include("solve_continuous.jl") end
@time @testset "Reachability.solve_hybrid_bouncing_ball" begin include("solve_hybrid_bouncing_ball.jl") end
@time @testset "Reachability.solve_hybrid_thermostat" begin include("solve_hybrid_thermostat.jl") end
