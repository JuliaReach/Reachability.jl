# Thermostat example
# See https://fenix.tecnico.ulisboa.pt/downloadFile/3779579688470/lygeros.pdf,
# Section 1.3.4.
# ============================

using HybridSystems, MathematicalSystems, LazySets, Polyhedra
import LazySets.HalfSpace
import Reachability.ReachSet

c_a = 0.1;
# Transition graph (automaton)
a = LightAutomaton(2);
add_transition!(a, 1, 2, 1);
add_transition!(a, 2, 1, 2);

# Mode on
A = hcat(-c_a);
B = hcat(30.);
X = HalfSpace([1.0], 22.0); # x <= 22
U = Singleton([c_a]);
m_on = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

# Mode off
A = hcat(-c_a);
B = hcat(0.0);
X = HalfSpace([-1.0], -18.0); # x >= 18
U = Singleton([0.0]);
m_off = ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U);

# Transition from on to off
A = hcat(1.0);
X = HPolytope([HalfSpace([-1.0], -21.0)]); # x >= 21
t_on2off = ConstrainedLinearDiscreteSystem(A, X)

# Transition from off to on
A = hcat(1.0);
X = HPolytope([HalfSpace([1.0], 19.0)]); # x <= 19
t_off2on = ConstrainedLinearDiscreteSystem(A, X)

m = [m_on, m_off];

r = [t_on2off, t_off2on];

# Switchings
s = [HybridSystems.AutonomousSwitching()];

HS = HybridSystem(a, m, r, s);

# initial condition in mode 1
X0 = Singleton([18.]);

system = InitialValueProblem(HS, X0);
input_options = Options(:mode=>"reach");
plot_vars = [0, 1]

problem_options = Options(:vars=>[1], :T=>5.0, :δ=>0.1, :plot_vars=>plot_vars,
                          :max_jumps=>1, :verbosity=>1,
                          :project_reachset => false);
options_input = merge(problem_options, input_options);
sol = solve(system, options_input);

# work-around for 1D plot
N = Float64
new_reach_sets = Vector{ReachSet{CartesianProductArray{N}, N}}(length(sol.Xk))
for (i, rs) in enumerate(sol.Xk)
    new_reach_sets[i] =
        ReachSet{CartesianProductArray{N}, N}(
            CartesianProductArray{N, HPolytope{N}}([rs.X]),
            rs.t_start, rs.t_end)
end
sol_processed = Reachability.ReachSolution(new_reach_sets, sol.options);
sol_proj = Reachability.ReachSolution(Reachability.project_reach(
    sol_processed.Xk, plot_vars, 1, sol.options), sol.options);
