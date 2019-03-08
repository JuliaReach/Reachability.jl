let # preventing global scope
    # ===================================================================
    # Discretization of a continuous-time system without input (ZeroSet)
    # ===================================================================
    A = sparse([1, 1, 2, 3, 4], [1, 2, 2, 4, 3], [1., 2., 3., 4., 5.], 4, 4)
    X0 = BallInf(zeros(4), 0.1)
    cont_sys_homog = IVP(LCS(A), X0)
    δ = 0.01

    # no bloating, use default options
    discr_sys_homog = discretize(cont_sys_homog, δ, algorithm="nobloating")
    @test inputdim(discr_sys_homog) == 0 # there is no input set!

    # no bloating, use Pade approximation
    discr_sys_homog = discretize(cont_sys_homog, δ, algorithm="nobloating", exp_method="pade")  
    @test inputdim(discr_sys_homog) == 0

    # bloating, default options
    discr_sys_homog = discretize(cont_sys_homog, δ)
    @test inputdim(discr_sys_homog) == 0

    # bloating, first order
    discr_sys_homog = discretize(cont_sys_homog, δ, algorithm="firstorder")
    @test inputdim(discr_sys_homog) == 0

    # bloating, forward interpolation
    discr_sys_homog = discretize(cont_sys_homog, δ, algorithm="forward")
    @test inputdim(discr_sys_homog) == 0

    # bloating, backward interpolation
    discr_sys_homog = discretize(cont_sys_homog, δ, algorithm="backward")
    @test inputdim(discr_sys_homog) == 0

    # bloating, forward interpolation using lazy operations (default)
    discr_sys_homog = discretize(cont_sys_homog, δ, algorithm="forward", set_operations="lazy")
    @test discr_sys_homog.x0 isa ConvexHull
    @test inputdim(discr_sys_homog) == 0

    # using concrete zonotopes
    discr_sys_homog = discretize(cont_sys_homog, δ, set_operations="zonotope")
    @test discr_sys_homog.x0 isa Zonotope
    @test inputdim(discr_sys_homog) == 0

    # ===============================================================
    # Discretization of a continuous-time system with constant input
    # ===============================================================
    U = ConstantInput(Ball2(ones(4), 0.5))
    B = Matrix{Float64}(1.0I, 4, 4)
    cont_sys = IVP(CLCCS(A, B, nothing, U), X0)

    # no bloating
    discr_sys = discretize(cont_sys, δ, algorithm="nobloating")
    @test inputdim(discr_sys) == 4

    inputs = next_set(inputset(discr_sys))
    @test dim(inputs) == 4
    @test isa(inputs, LinearMap)
    @test isa(inputs.X, Ball2) && inputs.X.center == ones(4) && inputs.X.radius == 0.5

    # no bloating, use Pade approximation
    discr_sys = discretize(cont_sys, δ, algorithm="nobloating", exp_method="pade")
    @test inputdim(discr_sys) == 4

    # bloating, use scaling and squaring method
    discr_sys = discretize(cont_sys, δ, exp_method="base")
    @test inputdim(discr_sys) == 4

    inputs = next_set(inputset(discr_sys))
    @test dim(inputs) == 4
    @test isa(inputs, MinkowskiSum)

    # bloating, use Pade approximation
    discr_sys = discretize(cont_sys, δ, exp_method="pade")
    @test inputdim(discr_sys) == 4

    discr_sys = discretize(cont_sys, δ, algorithm="firstorder")
    @test inputdim(discr_sys) == 4

    # using concrete zonotopes
    # since a Ball2 cannot be converted to a Zonotope, we use a box instead
    U = ConstantInput(BallInf(ones(4), 0.5))
    cont_sys = IVP(CLCCS(A, B, nothing, U), X0)
    discr_sys = discretize(cont_sys, δ, set_operations="zonotope")
    @test discr_sys.x0 isa Zonotope
    @test inputdim(discr_sys) == 4

    # ===================================================================
    # Discretization of a continuous-time system with time-varying input
    # ===================================================================
    Ui = [Ball2(0.01*i*ones(4), i*0.2) for i in 1:3]
    cont_sys = IVP(CLCCS(A, Matrix(1.0I, 4, 4), nothing, VaryingInput(Ui)), X0)

    # no bloating
    discr_sys = discretize(cont_sys, δ, algorithm="nobloating")
    Ui_d = inputset(discr_sys)
    @test length(Ui_d) == 3
    for (i, inputs) in enumerate(discr_sys.s.U)
        if i == 1
            @test dim(inputs) == 4
            @test isa(inputs, LinearMap)
            @test isa(inputs.X, Ball2) && inputs.X.center == 0.01*ones(4) && inputs.X.radius == 0.2
        elseif i == 2
            @test dim(inputs) == 4
            @test isa(inputs, LinearMap)
            @test isa(inputs.X, Ball2) && inputs.X.center == 0.01*2*ones(4) && inputs.X.radius == 0.2*2
        else
            @test dim(inputs) == 4
            @test isa(inputs, LinearMap)
            @test isa(inputs.X, Ball2) && inputs.X.center == 0.01*3*ones(4) && inputs.X.radius == 0.2*3
        end
    end

    # no bloating, use Pade approximation
    discr_sys = discretize(cont_sys, δ, algorithm="nobloating", exp_method="pade")

    # bloating
    discr_sys = discretize(cont_sys, δ, exp_method="base")
    Ui_d = inputset(discr_sys)
    @test length(Ui_d) == 3
    for (i, inputs) in enumerate(discr_sys.s.U)
        if i == 1
            @test dim(inputs) == 4
            @test isa(inputs, MinkowskiSum)
            X = inputs.X.X
            @test isa(X, Ball2) && X.center == 0.01*ones(4) && X.radius == 0.2
        elseif i == 2
            @test dim(inputs) == 4
            @test isa(inputs, MinkowskiSum)
            X = inputs.X.X
            @test isa(X, Ball2) && X.center == 0.01*2*ones(4) && X.radius == 0.2*2
        else
            @test dim(inputs) == 4
            @test isa(inputs, MinkowskiSum)
            X = inputs.X.X
            @test isa(X, Ball2) && X.center == 0.01*3*ones(4) && X.radius == 0.2*3
        end
    end

    # bloating, use Pade approximation
    discr_sys = discretize(cont_sys, δ, exp_method="pade")

    # lazy symmetric interval hull
    discr_sys = discretize(cont_sys, δ, sih_method="lazy")

    # concrete symmetric interval hull
    discr_sys = discretize(cont_sys, δ, sih_method="concrete")
end
