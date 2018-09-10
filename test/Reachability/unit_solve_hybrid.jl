# Bouncing ball example
# See: https://juliareach.github.io/SX.jl/latest/examples/bball.html
# ============================

using HybridSystems, MathematicalSystems, LazySets

# Transition graph (automaton)
a = LightAutomaton(1)
add_transition!(a, 1, 1, 1);

# Modes
A = [0.0 1.0; 0.0 0.0]
B = reshape([0.0, -1.0], (2, 1))
X = [HalfSpace([1.0, 0.0], 0.0)]
U = Singleton([1.0])
m = [ConstrainedLinearControlContinuousSystem(A, eye(size(B, 1)), X, B*U)];

# Reset maps
A = [1.0 0.0; 0.0 -0.75]
B = zeros(2, 0)

X = [HalfSpace([0.0, 1.0], 0.0),
           HalfSpace([-1.0, 0.0], 0.0),# x = 0
           HalfSpace([1.0, 0.0], 0.0)] # v < 0
U = Vector{LazySet{Float64}}()
r = [ConstrainedLinearControlDiscreteSystem(A, B, X, U)];

# Switchings
s = [HybridSystems.AutonomousSwitching()];

H = HybridSystem(a, m, r, s)

# initial condition in mode 1
X0 = BallInf(zeros(2), 0.01)

# calculate reachable states up to time T
prob = InitialValueProblem(H, X0)
input_options = Options(:mode=>"reach")


problem_options = Options(:vars=>[1,2],:T=>1.0,
                    :δ=>0.01, :plot_vars => [1, 2]);
options = merge(problem_options, input_options)
sol = solve_hybrid(H, X0, options)
