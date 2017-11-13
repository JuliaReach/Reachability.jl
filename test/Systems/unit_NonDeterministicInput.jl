import Reachability.Systems: ConstantNonDeterministicInput, TimeVaryingNonDeterministicInput

# Testing null constant input
state = start(ConstantNonDeterministicInput())
@test isa(state.set, VoidSet)

# Testing constant non-deterministic input
c = zeros(4); r = 0.1
U = BallInf(c, r)
input_iter =  ConstantNonDeterministicInput(U)

@test length(input_iter) == 1
state = start(input_iter)
input_remains_constant = state.set.center == c && state.set.radius == r
@test input_remains_constant
for i in 1:5
    state = next(input_iter, state)
    input_remains_constant = input_remains_constant && state.set.center == c && state.set.radius == r
end

# Testing time-varying non-deterministic input
c = [zeros(4), 1. + ones(4), 2. + ones(4)]
r = [0.1, 0.2, 0.3]

input_iter =  TimeVaryingNonDeterministicInput([BallInf(c[i], r[i]) for i in 1:3])

@test length(input_iter) == 3
state = start(input_iter)
input_is_varying = state.set.center == c[1] && state.set.radius == r[1]
@test input_is_varying
for i in 2:3
    state = next(input_iter, state)
    input_is_varying = input_is_varying && state.set.center == c[i] && state.set.radius == r[i]
end
@test input_is_varying
