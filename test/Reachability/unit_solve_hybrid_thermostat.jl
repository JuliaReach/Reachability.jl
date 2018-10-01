# Bouncing ball example
# See: https://juliareach.github.io/SX.jl/latest/examples/bball.html
# ============================

using HybridSystems, MathematicalSystems, LazySets, Plots, Polyhedra

c_a = 0.1;
# Transition graph (automaton)
a = LightAutomaton(2);
add_transition!(a, 1, 2, 1);
add_transition!(a, 2, 1, 2);

# Mode off
A = hcat(-c_a);
B = hcat(0.0);
X = HPolytope([HalfSpace([1.0], 22.0)]); # x <= 22
U = Singleton([0.0]);
m_on = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

# Mode on
A = hcat(-c_a);
B = hcat(30.);
X = HPolytope([HalfSpace([-1.0], -18.0)]); # x >= 18
U = Singleton([c_a]);
m_off = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

# Transition from on to off
A_off = hcat(1.0);
X_off = HPolytope([HalfSpace([-1.0], -21.0)]); # x >= 21

# Transition from off to on
A_on = hcat(1.0);
X_on = HPolytope([HalfSpace([1.0], 19.0)]); # x <= 19

m = [m_on, m_off];

r = [ConstrainedLinearDiscreteSystem(A_on, X_on),
     ConstrainedLinearDiscreteSystem(A_off, X_off)];

# Switchings
s = [HybridSystems.AutonomousSwitching()];

HS = HybridSystem(a, m, r, s);

# initial condition in mode 1
X0 = Singleton([18.]);

# calculate reachable states up to time T
prob = InitialValueProblem(HS, X0);
input_options = Options(:mode=>"reach");
plot_vars = [0, 1]

problem_options = Options(:vars=>[1], :T=>10.0, :δ=>0.01, :verbosity=>1, :plot_vars=>plot_vars);
options_input = merge(problem_options, input_options);
sol = solve_hybrid(HS, X0, options_input);

# work-around for 1D plot
N = Float64
sol_processed = ReachSolution([CartesianProductArray{N, HPolytope{N}}([x]) for x in sol.Xk], sol.options);
sol_proj = Reachability.project_reach(plot_vars, 1, options_input.dict[:δ], sol_processed.Xk, "");

plot(sol_proj)
