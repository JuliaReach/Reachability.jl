# Bouncing ball example
# See: https://juliareach.github.io/SX.jl/latest/examples/bball.html
# ============================

using HybridSystems, MathematicalSystems, LazySets, Plots

# Transition graph (automaton)
a = LightAutomaton(2)
add_transition!(a, 1, 2, 1);
add_transition!(a, 2, 1, 2);

# Mode on
A = [-0.1]
B = [30]
X = HPolytope([HalfSpace([1.0, 0.0], 22.0)]) # x <= 22
U = Singleton([0.1])
m = [ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U)];

# Mode on
A = [-0.1]
B = [30]
X = HPolytope([HalfSpace([-1.0, 0.0], 18.0)]) # x >= 18
U = Singleton([0.1])
m = [ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U)];

# Transition from on to off
A_off = [1.0 0.0]
X_off = HPolytope([HalfSpace([-21.0, 21.0], 0.0)) # x >= 21

# Transition from off to on
A_on = [1.0 0.0]
X_on = HPolytope([HalfSpace([0.0, 19.0], 0.0),  # x <= 19



r = [ConstrainedLinearDiscreteSystem(A_on, X_on), ConstrainedLinearDiscreteSystem(A_off, X_off)];

# Switchings
s = [HybridSystems.AutonomousSwitching()];

HS = HybridSystem(a, m, r, s)

# initial condition in mode 1
X0 = Singleton([18])

# calculate reachable states up to time T
prob = InitialValueProblem(HS, X0)
input_options = Options(:mode=>"reach")

problem_options = Options(:vars=>[1], :T=>10.0, :δ=>0.01, :verbosity=>1);
options_input = merge(problem_options, input_options)
using Polyhedra
sol = solve_hybrid(HS, X0, options_input);

plot(sol, indices=1:2:length(sol.Xk))
