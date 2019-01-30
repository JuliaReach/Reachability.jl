# Bouncing ball example
# See: https://juliareach.github.io/SX.jl/latest/examples/bball.html
# ============================

using Reachability, HybridSystems, MathematicalSystems, LazySets, Polyhedra
import LazySets.HalfSpace
import Reachability.constrained_dimensions

# Transition graph (automaton)
a = LightAutomaton(1);
add_transition!(a, 1, 1, 1);

# mode
A = [0.0 1.0; 0.0 0.0];
B = reshape([0.0, -1.0], (2, 1));
X = HalfSpace([-1.0, 0.0], 0.0); # x >= 0
U = Singleton([1.0]);
m = [ConstrainedLinearControlContinuousSystem(A, Matrix{eltype(A)}(I, size(B, 1), size(B, 1)), X, B*U)];

# reset map
A = [1.0 0.0; 0.0 -0.75];
X = HPolyhedron([HalfSpace([0.0, 1.0], 0.0),   # v <= 0
                 HalfSpace([-1.0, 0.0], 0.0),  # x >= 0
                 HalfSpace([1.0, 0.0], 0.0)]); # x <= 0
r = [ConstrainedLinearDiscreteSystem(A, X)];

# switching
s = [HybridSystems.AutonomousSwitching()];

HS = HybridSystem(a, m, r, s);

# initial condition in mode 1
X0 = Hyperrectangle(low=[10, 0.0], high=[10.2, 0.0]);

inits = [(1, X0)]

system = InitialValueProblem(HS, inits);

options = Options(:mode=>"reach", :vars=>[1, 2], :T=>5.0, :δ=>0.1,
                  :plot_vars=>[1, 2], :max_jumps=>1, :verbosity=>1);


# test constrained_dimensions
@test constrained_dimensions(HS) == Dict{Int,Vector{Int}}(1=>[1, 2])

# default algorithm
sol = solve(system, options);

# specify lazy discrete post operator
sol = solve(system, options, BFFPSV18(), LazyDiscretePost());

# specify overapproximating discrete post operator
sol = solve(system, options, BFFPSV18(), ApproximatingDiscretePost());
