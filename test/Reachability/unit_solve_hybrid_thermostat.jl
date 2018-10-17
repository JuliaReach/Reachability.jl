# Thermostat example
# See https://fenix.tecnico.ulisboa.pt/downloadFile/3779579688470/lygeros.pdf,
# Section 1.3.4.
# ============================

using Reachability, HybridSystems, MathematicalSystems, LazySets, Polyhedra, Optim
import LazySets.HalfSpace
import Reachability.ReachSet

c_a = 0.1;

# transition graph (automaton)
a = LightAutomaton(2);
add_transition!(a, 1, 2, 1);
add_transition!(a, 2, 1, 2);

# mode on
A = hcat(-c_a);
B = hcat(30.);
X = HalfSpace([1.0], 22.0); # x <= 22
U = Singleton([c_a]);
m_on = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

# mode off
A = hcat(-c_a);
B = hcat(0.0);
X = HalfSpace([-1.0], -18.0); # x >= 18
U = Singleton([0.0]);
m_off = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

# transition from on to off
A = hcat(1.0);
X = HalfSpace([-1.0], -21.0); # x >= 21
t_on2off = ConstrainedLinearDiscreteSystem(A, X)

# transition from off to on
A = hcat(1.0);
X = HalfSpace([1.0], 19.0); # x <= 19
t_off2on = ConstrainedLinearDiscreteSystem(A, X)

m = [m_on, m_off];

r = [t_on2off, t_off2on];

# switchings
s = [HybridSystems.AutonomousSwitching()];

HS = HybridSystem(a, m, r, s);

# initial condition in mode off
X0 = Singleton([18.]);
inits = [(2, X0)]

system = InitialValueProblem(HS, X0);

options = Options(:mode=>"reach", :vars=>[1], :T=>5.0, :δ=>0.1,
                  :max_jumps=>1, :verbosity=>1);

# default algorithm
sol = solve(system, options);

# specify lazy discrete post operator
sol = solve(system, options, Reachability.BFFPSV18(),
            Reachability.ReachSets.LazyTextbookDiscretePost());

# specify overapproximating discrete post-operator algorithm
sol = solve(system, options, Reachability.BFFPSV18(),
            Reachability.ReachSets.ApproximatingDiscretePost());
