# Bouncing ball example
# See: https://juliareach.github.io/SX.jl/latest/examples/bball.html
# ============================

using HybridSystems, MathematicalSystems, LazySets, Polyhedra
import LazySets.HalfSpace

# Transition graph (automaton)
a = LightAutomaton(1);
add_transition!(a, 1, 1, 1);

# Modes
A = [0.0 1.0; 0.0 0.0];
B = reshape([0.0, -1.0], (2, 1));
X = HalfSpace([-1.0, 0.0], 0.0); # x >= 0
U = Singleton([1.0]);
m = [ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U)];

# Reset maps
A = [1.0 0.0; 0.0 -0.75];
X = HPolytope([HalfSpace([0.0, 1.0], 0.0),   # v <= 0
               HalfSpace([-1.0, 0.0], 0.0),  # x >= 0
               HalfSpace([1.0, 0.0], 0.0)]); # x <= 0
r = [ConstrainedLinearDiscreteSystem(A, X)];

# Switchings
s = [HybridSystems.AutonomousSwitching()];

HS = HybridSystem(a, m, r, s);

# initial condition in mode 1
X0 = Hyperrectangle(low=[10, 0.0], high=[10.2, 0.0]);

prob = InitialValueProblem(HS, X0);
input_options = Options(:mode=>"reach");

problem_options = Options(:vars=>[1,2], :T=>5.0, :δ=>0.1, :plot_vars=>[1, 2],
                          :max_jumps=>1, :verbosity=>1);
options_input = merge(problem_options, input_options);
sol = solve(HS, X0, options_input);
