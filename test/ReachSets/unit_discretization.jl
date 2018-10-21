import Reachability.ReachSets.discretize

let # preventing global scope
    # ===================================================================
    # Discretization of a continuous-time system without input (ZeroSet)
    # ===================================================================
    A = sparse([1, 1, 2, 3, 4], [1, 2, 2, 4, 3], [1., 2., 3., 4., 5.], 4, 4)
    X0 = BallInf(zeros(4), 0.1)
    #U = ZeroSet(size(A, 1))
    cont_sys_homog = ContinuousSystem(A, X0)
    δ = 0.01

    # no bloating, do not use Pade approximation
    discr_sys_homog = discretize(cont_sys_homog, δ, approx_model="nobloating", pade_expm=false)
    @test inputdim(discr_sys_homog) == 0 # there is no input set!

    # no bloating, use Pade approximation
    discr_sys_homog = discretize(cont_sys_homog, δ, approx_model="nobloating", pade_expm=true)
    @test inputdim(discr_sys_homog) == 0

    # bloating, do not use Pade approximation
    discr_sys_homog = discretize(cont_sys_homog, δ, pade_expm=false)
    @test inputdim(discr_sys_homog) == 0

    # bloating, first order
    discr_sys_homog = discretize(cont_sys_homog, δ, approx_model="firstorder")
    @test inputdim(discr_sys_homog) == 0

    # ===============================================================
    # Discretization of a continuous-time system with constant input
    # ===============================================================
    U = Ball2(ones(4), 0.5)
    cont_sys = ContinuousSystem(A, X0, U)

    # no bloating, do not use Pade approximation
    discr_sys = discretize(cont_sys, δ, approx_model="nobloating", pade_expm=false)
    @test inputdim(discr_sys) == 4
    #@test length(discr_sys.U) == 1 # a constant input has no length method!
    inputs = next_set(inputset(discr_sys)) # getter: inputset(discr_sys) = discr_sys.s.U
    @test dim(inputs) == 4
    @test isa(inputs, LinearMap)
    @test isa(inputs.X, Ball2) && inputs.X.center == ones(4) && inputs.X.radius == 0.5

    # no bloating, use Pade approximation
    discr_sys = discretize(cont_sys, δ, approx_model="nobloating", pade_expm=true)
    @test inputdim(discr_sys) == 4

    # bloating, do not use Pade approximation
    discr_sys = discretize(cont_sys, δ, pade_expm=false)
    @test inputdim(discr_sys) == 4

    #@test length(discr_sys.U) == 1
    inputs = next_set(inputset(discr_sys))
    @test dim(inputs) == 4
    @test isa(inputs, MinkowskiSum)

    # bloating, use Pade approximation
    discr_sys = discretize(cont_sys, δ, pade_expm=true)
    @test inputdim(discr_sys) == 4

    discr_sys = discretize(cont_sys, δ, approx_model="firstorder")
    @test inputdim(discr_sys) == 4

    # ===================================================================
    # Discretization of a continuous-time system with time-varying input
    # ===================================================================
    Ui = [Ball2(0.01*i*ones(4), i*0.2) for i in 1:3]
    cont_sys = ContinuousSystem(A, X0, Ui)

    # no bloating, do not use Pade approximation
    discr_sys = discretize(cont_sys, δ, approx_model="nobloating", pade_expm=false)
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
    discr_sys = discretize(cont_sys, δ, approx_model="nobloating", pade_expm=true)

    # bloating, do not use Pade approximation
    discr_sys = discretize(cont_sys, δ, pade_expm=false)
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
    discr_sys = discretize(cont_sys, δ, pade_expm=true)

    # lazy symmetric interval hull
    discr_sys = discretize(cont_sys, δ, lazy_sih=true)
    discr_sys = discretize(cont_sys, δ, lazy_sih=false)
end
