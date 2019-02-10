# =================================================================
# Bouncing-ball model
# See https://juliareach.github.io/SX.jl/latest/examples/bball.html
# =================================================================

using Reachability, HybridSystems, MathematicalSystems, LazySets, LinearAlgebra
using LazySets: HalfSpace  # resolve name-space conflicts with Polyhedra

function bouncing_ball()
    # automaton structure
    automaton = LightAutomaton(1)

    # mode 1
    A = [0.0 1.0; 0.0 0.0]
    B = reshape([0.0, -1.0], (2, 1))
    U = Singleton([1.0])
    inv = HalfSpace([-1.0, 0.0], 0.0) # x >= 0
    m1 = ConstrainedLinearControlContinuousSystem(
            A, Matrix{eltype(A)}(I, size(B, 1), size(B, 1)), inv, B*U)

    # modes
    modes = [m1]

    # transition from m1 to m1 (self-loop)
    add_transition!(automaton, 1, 1, 1)
    A = [1.0 0.0; 0.0 -0.75]
    guard = HPolyhedron([HalfSpace([0.0, 1.0], 0.0),   # v ≤ 0
                         HalfSpace([-1.0, 0.0], 0.0),  # x ≥ 0
                         HalfSpace([1.0, 0.0], 0.0)])  # x ≤ 0
    t1 = ConstrainedLinearMap(A, guard)

    # transition annotations
    resetmaps = [t1]

    # switching
    switchings = [AutonomousSwitching()]

    ℋ = HybridSystem(automaton, modes, resetmaps, switchings)

    # initial condition in mode 1
    X0 = Hyperrectangle(low=[10.0, 0.0], high=[10.2, 0.0])
    initial_condition = [(1, X0)]

    system = InitialValueProblem(ℋ, initial_condition)

    options = Options(:mode=>"reach", :vars=>[1, 2], :T=>5.0, :δ=>0.1,
                      :plot_vars=>[1, 2], :max_jumps=>1, :verbosity=>1)

    return (system, options)
end
