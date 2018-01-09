import Reachability.Systems: ConstantNonDeterministicInput, TimeVaryingNonDeterministicInput

# Testing constant non-deterministic input
c = zeros(4); r = 0.1
U = BallInf(c, r)
U = ConstantNonDeterministicInput(U)

@test length(U) == 1
inputs = next_set(U)
input_remains_constant = inputs.center == c && inputs.radius == r
@test input_remains_constant
for i in 1:5
    inputs = next(U, i)[1]
    input_remains_constant = input_remains_constant && inputs.center == c && inputs.radius == r
end

# Testing time-varying non-deterministic input
c = [zeros(4), 1. + ones(4), 2. + ones(4)]
r = [0.1, 0.2, 0.3]

U = TimeVaryingNonDeterministicInput([BallInf(c[i], r[i]) for i in 1:3])

@test length(U) == 3
input_is_varying = true
i = 1
for u in U
    input_is_varying = input_is_varying && u.center == c[i] && u.radius == r[i]
    i += 1
end
@test input_is_varying
