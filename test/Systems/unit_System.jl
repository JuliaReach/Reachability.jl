# ======================================================
# Testing continuous-time homogeneous system (no input)
# ======================================================
A = sparse([1, 1, 2, 3, 4], [1, 2, 2, 4, 3], [1., 2., 3., 4., 5.], 4, 4)
X0 = BallInf(zeros(4), 0.1)
cont_sys_homog = ContinuousSystem(A, X0)

# Check if the input is constant
@test isa(cont_sys_homog.U, Systems.ConstantNonDeterministicInput)
# Check if the input is empty
@test isa(start(cont_sys_homog.U).sf, VoidSet)
# Check data fields
@test cont_sys_homog.A == A
@test cont_sys_homog.X0.center == zeros(4) && cont_sys_homog.X0.radius == 0.1

# ===================================================
# Testing continuous-time system with constant input
# ===================================================
U = Ball2(ones(4), 0.5)
cont_sys = ContinuousSystem(A, X0, U)

# check initial state
@test cont_sys.X0.center ≈ zeros(4) && cont_sys.X0.radius ≈ 0.1

# recover input
input_state = start(cont_sys.U)

@test input_state.sf.center ≈ ones(4) && input_state.sf.radius ≈ 0.5

# ========================================================
# Testing continuous-time system with time-varying input
# ========================================================
Ui = [Ball2(0.01*i*ones(4), i*0.2) for i in 1:3]
cont_sys = ContinuousSystem(A, X0, Ui)

input_state = start(cont_sys.U)
@test input_state.sf.center ≈ 0.01*ones(4) && input_state.sf.radius ≈ 0.2

input_state = next(cont_sys.U, input_state)
@test input_state.sf.center ≈ 0.02*ones(4) && input_state.sf.radius ≈ 0.4

input_state = next(cont_sys.U, input_state)
@test input_state.sf.center ≈ 0.03*ones(4) && input_state.sf.radius ≈ 0.6

# =========================================
# Testing discrete-time homogeneous system
# =========================================
δ = 0.01
discr_sys_homog = DiscreteSystem(A, X0, δ)

# Check if the input is constant
@test isa(discr_sys_homog.U, Systems.ConstantNonDeterministicInput)
# Check if the input is empty
@test isa(start(discr_sys_homog.U).sf, VoidSet)
# Check data fields
@test discr_sys_homog.A == A
@test discr_sys_homog.X0.center == zeros(4) && discr_sys_homog.X0.radius == 0.1
@test discr_sys_homog.δ == δ

# =================================================
# Testing discrete-time system with constant input
# =================================================
U = Ball2(ones(4), 0.5)
discr_sys = DiscreteSystem(A, X0, δ, U)

# check initial state
@test discr_sys.X0.center ≈ zeros(4) && discr_sys.X0.radius ≈ 0.1

# recover input
input_state = start(discr_sys.U)

@test input_state.sf.center ≈ ones(4) && input_state.sf.radius ≈ 0.5

# =====================================================
# Testing discrete-time system with time-varying input
# =====================================================
Ui = [Ball2(0.01*i*ones(4), i*0.2) for i in 1:3]
discr_sys = DiscreteSystem(A, X0, δ, Ui)

input_state = start(discr_sys.U)
@test input_state.sf.center ≈ 0.01*ones(4) && input_state.sf.radius ≈ 0.2

input_state = next(discr_sys.U, input_state)
@test input_state.sf.center ≈ 0.02*ones(4) && input_state.sf.radius ≈ 0.4

input_state = next(discr_sys.U, input_state)
@test input_state.sf.center ≈ 0.03*ones(4) && input_state.sf.radius ≈ 0.6
