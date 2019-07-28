# ==============================================================================
# Filtered-oscillator model
# See: https://flowstar.org/benchmarks/filtered-oscillator/
#
# Modification: made parametric in variables x and y so they can be changed.
# We require that y > x for simpler computation of `filter1`.
# ==============================================================================

using HybridSystems, MathematicalSystems, LazySets, LinearAlgebra, Reachability
import LazySets.HalfSpace  # resolve name-space conflicts with Polyhedra

function filtered_oscillator(n0::Int=6, x=1, y=2)
    n = n0 + 2

    # transition graph (automaton)
    a = LightAutomaton(4)
    add_transition!(a, 3, 4, 1)
    add_transition!(a, 4, 2, 2)
    add_transition!(a, 2, 1, 3)
    add_transition!(a, 1, 3, 4)

    # common flow
    A = zeros(n, n)
    for i = 2:n-1
        A[i,i-1] = 5.
        A[i,i] = -5.
    end
    A[x, x] = -2.
    A[y, y] = -1.
    if x > 1
        filter1 = 1
    elseif y > 2
        filter1 = 2
    else
        filter1 = 3
    end
    A[filter1,x], A[filter1, filter1] = 5., -5.

    # commond constraints
    hs1_a = zeros(n)
    hs1_a[x] = -0.714286
    hs1_a[y] = -1.0
    hs1 = HalfSpace(hs1_a, 0.0)  # 0.714286*x + y >= 0
    hs2_a = zeros(n)
    hs2_a[x] = 1.0
    hs2 = HalfSpace(hs2_a, 0.0)  # x <= 0
    hs3_a = zeros(n)
    hs3_a[x] = 0.714286
    hs3_a[y] = 1.0
    hs3 = HalfSpace(hs3_a, 0.0)  # 0.714286*x + y <= 0
    hs4_a = zeros(n)
    hs4_a[x] = -1.0
    hs4 = HalfSpace(hs4_a, 0.0)  # x >= 0

    # modes

    # mode 1
    b = zeros(n)
    b[x] = 1.4
    b[y] = -0.7
    X = HPolyhedron([hs1, hs2])
    m_1 = CACS(A, b, X)

    # mode 2
    b = zeros(n)
    b[x] = -1.4
    b[y] = 0.7
    X = HPolyhedron([hs2, hs3])
    m_2 = CACS(A, b, X)

    # mode 3
    b = zeros(n)
    b[x] = 1.4
    b[y] = -0.7
    X = HPolyhedron([hs4, hs1])
    m_3 = CACS(A, b, X)

    # mode 4
    b = zeros(n)
    b[x] = -1.4
    b[y] = 0.7
    X = HPolyhedron([hs3, hs4])
    m_4 = CACS(A, b, X)

    m = [m_1, m_2, m_3, m_4]

    # transitions

    # transition l3 -> l4
    X_l3l4 = HPolyhedron([hs4, hs1, hs3])
    r1 = ConstrainedIdentityMap(n, X_l3l4)

    # transition l4 -> l2
    X_l4l2 = HPolyhedron([hs3, hs4, hs2])
    r2 = ConstrainedIdentityMap(n, X_l4l2)

    # transition l2 -> l1
    X_l2l1 = HPolyhedron([hs2, hs1, hs3])
    r3 = ConstrainedIdentityMap(n, X_l2l1)

    # transition l1 -> l3
    X_l1l3 = HPolyhedron([hs1, hs4, hs2])
    r4 = ConstrainedIdentityMap(n, X_l1l3)

    r = [r1, r2, r3, r4]

    # switchings
    s = [HybridSystems.AutonomousSwitching()]

    HS = HybridSystem(a, m, r, s)

    # initial condition in mode 3
    low = zeros(n)
    low[x] = 0.2
    low[y] = -0.1
    high = zeros(n)
    high[x] = 0.3
    high[y] = 0.1
    X0 = Hyperrectangle(low=low, high=high)

    problem = InitialValueProblem(HS, [(3, X0)])

    options = Options(:T=>Inf, :mode=>"reach", :max_jumps=>5, :verbosity=>1)

    return (problem, options)
end
