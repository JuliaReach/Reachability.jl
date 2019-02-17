using Reachability.ReachSets: BFFPSV18

#===== Projectile model =====
We test the line plot!(x->x, x->-24*x+375, 0., 50., line=2., color="red", linestyle=:solid, legend=:none)
=#
A = [0. 0.5 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.7 ; 0. 0. 0. 0.]
X0 = Singleton([0.,5.,100.,0])
U = Singleton([0.,0.,0.,-9.81])
S = ContinuousSystem(A, X0, U)
time_horizon = 20.0
prec = 1e-4
initial_δ = 0.5
algorithm(N, δ) =
    solve(S,
          Options(:mode => "check", :plot_vars => [1, 3], :T => time_horizon,
                  :property=>LinearConstraintProperty([24., 0., 1, 0],  375.));
          op=BFFPSV18(:δ => δ, :vars => [1, 3], :partition=>[1:2, 3:4])
         ).satisfied

Reachability.tune_δ(algorithm, time_horizon, prec, initial_δ)
