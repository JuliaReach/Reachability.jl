# Simple example of 2 elementary dynamics to test low-dimensional discrete post
# ============================

using Reachability, HybridSystems, MathematicalSystems, LazySets, Polyhedra

# fix namamespace conflicts with Polyhedra
import LazySets.HalfSpace
import Reachability.ReachSet
import Reachability.constrained_dimensions

n = 2
system_dimension = n + 2
z = zeros(n)
# transition graph (automaton)
a = LightAutomaton(2);
add_transition!(a, 1, 2, 1);
add_transition!(a, 2, 1, 2);

# common B and U
U = Singleton([1.0]);
B = [0.0; 0.0; z];

# mode1
A = zeros(system_dimension, system_dimension)
A[1,1], A[2,1] = 1.0 , 1.0
X = HalfSpace([0.0; 1.0; z], 3.0); # y <= 3
r_1 = ConstrainedLinearControlContinuousSystem(A, Matrix{eltype(A)}(I, size(B, 1), size(B, 1)), X, B*U);

# mode2
A = zeros(system_dimension, system_dimension)
A[1,1], A[2,2] = 1.0 , 1.0
X = HalfSpace([0.0; 1.0; z], 10.0); # x <= 10
r_2 = ConstrainedLinearControlContinuousSystem(A, Matrix{eltype(A)}(I, size(B, 1), size(B, 1)), X, B*U);

# common resets
A_trans = eye(system_dimension)
# transition from on to off
X = HalfSpace([-1.0; 0.0; z], -8.0); # x >= 8
t_r1tor2 = ConstrainedLinearDiscreteSystem(A, X)

# transition from on to off
X = HalfSpace([-1.0; 0.0; z], -11.0); # x >= 11
t_r2tor1 = ConstrainedLinearDiscreteSystem(A, X)

m = [r_1, r_2];

r = [t_r1tor2,t_r2tor1];

# switchings
s = [HybridSystems.AutonomousSwitching()];

HS = HybridSystem(a, m, r, s);

# initial condition in mode к1
X0 = Singleton([0.0; 0.0; z]);
inits = [(1, X0)]

system = InitialValueProblem(HS, X0);

options = Options(:mode=>"reach", :vars=>1:system_dimension, :T=>5.0, :δ=>0.1,
                  :plot_vars=>[1, 2], :max_jumps=>1, :verbosity=>1);

# default algorithm
sol = solve(system, options);

# specify lazy discrete post operator
sol = solve(system, options, BFFPSV18(), LazyDiscretePost());

# specify overapproximating discrete post-operator algorithm
sol = solve(system, options, BFFPSV18(), ApproximatingDiscretePost());
