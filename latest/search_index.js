var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Reachability.jl-1",
    "page": "Home",
    "title": "Reachability.jl",
    "category": "section",
    "text": "Reachability is a Julia package for approximating the reachable states and checking safety properties of affine systems."
},

{
    "location": "index.html#Reachable-states-approximation-1",
    "page": "Home",
    "title": "Reachable states approximation",
    "category": "section",
    "text": "In a nutshell, we overapproximate the reachable states of an affine system by solving a set-based recurrence. The key idea is that we first decompose the system into (low-dimensional) subsystems and later compose the results as a Cartesian product. Thus we have to solve many cheap problems instead of one hard problem. Since solving the recurrence scales superlinearly with the dimension, this approach is very scalable."
},

{
    "location": "index.html#Decomposition-error-1",
    "page": "Home",
    "title": "Decomposition error",
    "category": "section",
    "text": "Decomposition typically involves a loss in precision, and so does this approach. The good thing is that we can decompose the recurrence as well, which allows us to analyze each of the subsystems independently by only referring to the initial states of the other subsystems. Consequently, there are two main sources for precision loss:Decomposition of the initial states: If two subsystems are interdependent initially.\nRepresentation of the reachable states as a Cartesian product: If two subsystems are interdependent in the dynamics.\nRepresentation of the reachable states in general: The reachable states of affine systems cannot be represented precisely in all cases. This is a problem that all approaches suffer from. We overapproximate the reachable states by (unions of) convex polytopes."
},

{
    "location": "index.html#Checking-safety-properties-1",
    "page": "Home",
    "title": "Checking safety properties",
    "category": "section",
    "text": "The problem of checking a safety property can be reduced to a reachability problem. We provide special support for this reduction by inlining the property check into the reachable states computation. This has two benefits:We fail fast when the property is violated in our abstraction.\nThe check is usually cheaper than computing the full reachable states. This is because we are often only interested in an upper or lower bound of a variable."
},

{
    "location": "index.html#Lazy-sets-1",
    "page": "Home",
    "title": "Lazy sets",
    "category": "section",
    "text": "To represent sets of states, we use the LazySets package which provides exact but lazy (i.e. symbolic) representations of common sets."
},

{
    "location": "index.html#Library-outline-1",
    "page": "Home",
    "title": "Library outline",
    "category": "section",
    "text": "Pages = [\n    \"lib/interface.md\",\n    \"lib/systems.md\",\n    \"lib/transformations.md\",\n    \"lib/discretize.md\",\n    \"lib/distributed.md\"\n]\nDepth = 2"
},

{
    "location": "lib/interface.html#",
    "page": "User interface",
    "title": "User interface",
    "category": "page",
    "text": ""
},

{
    "location": "lib/interface.html#User-interface-1",
    "page": "User interface",
    "title": "User interface",
    "category": "section",
    "text": "This section of the manual describes the main user interface.Pages = [\"interface.md\"]CurrentModule = Reachability"
},

{
    "location": "lib/interface.html#Reachability.solve",
    "page": "User interface",
    "title": "Reachability.solve",
    "category": "function",
    "text": "solve(system, options)  or  solve(system, :key1 => val1, [...], keyK => valK)\n\nSolves a reachability problem s.t. the given options. If some options are not defined, we may fall back to default values.\n\nInput\n\nsystem  – a (discrete or continuoues) system specification\noptions – options for solving the problem\n\nOutput\n\nA solution object whose content depends on the input options.\n\nNotes\n\nTo see all available input options, see keys(Reachability.available_keywords.dict).\n\n\n\n"
},

{
    "location": "lib/interface.html#Posing-and-solving-a-reachability-problem-1",
    "page": "User interface",
    "title": "Posing and solving a reachability problem",
    "category": "section",
    "text": "A reachability problem is characterized by an AbstractSystem together with an Options structure.solve"
},

{
    "location": "lib/systems.html#",
    "page": "Systems",
    "title": "Systems",
    "category": "page",
    "text": ""
},

{
    "location": "lib/systems.html#Systems-1",
    "page": "Systems",
    "title": "Systems",
    "category": "section",
    "text": "This module provides types to represent systems of affine ODEs in both discrete and continuous time.Pages = [\"systems.md\"]\nDepth = 3CurrentModule = Reachability.Systems"
},

{
    "location": "lib/systems.html#Types-of-systems-1",
    "page": "Systems",
    "title": "Types of systems",
    "category": "section",
    "text": "MathematicalSystems.jl provides some convenience types and methods to work with mathematical systems models. Every system inherits from AbstractSystem.We support the following two concrete types of systems."
},

{
    "location": "lib/systems.html#Discrete-system-1",
    "page": "Systems",
    "title": "Discrete system",
    "category": "section",
    "text": "A discrete system consists of a matrix representing the system dynamics, a set of initial states, a set of nondeterministic inputs, and a discretization step δ."
},

{
    "location": "lib/systems.html#Continuous-system-1",
    "page": "Systems",
    "title": "Continuous system",
    "category": "section",
    "text": "A continuous system consists of a matrix representing the system dynamics, a set of initial states, and a set of nondeterministic inputs."
},

{
    "location": "lib/systems.html#Nondeterministic-inputs-1",
    "page": "Systems",
    "title": "Nondeterministic inputs",
    "category": "section",
    "text": "The above systems may contain nondeterministic inputs, which are wrapped in special types. Every nondeterministic input representation inherits from NonDeterministicInput.The inputs are closely related to a DiscreteSystem in the sense that for each discrete time step the input set may change. We support iteration through the inputs over time."
},

{
    "location": "lib/systems.html#Constant-nondeterministic-inputs-1",
    "page": "Systems",
    "title": "Constant nondeterministic inputs",
    "category": "section",
    "text": "Constant nondeterministic inputs are chosen from a set of values that does not change over time. Note that, while the set is constant, the inputs themselves vary over time."
},

{
    "location": "lib/systems.html#Time-varying-nondeterministic-inputs-1",
    "page": "Systems",
    "title": "Time-varying nondeterministic inputs",
    "category": "section",
    "text": "Time-varying nondeterministic inputs are chosen from a set of values that changes over time (with each time step)."
},

{
    "location": "lib/transformations.html#",
    "page": "Transformations",
    "title": "Transformations",
    "category": "page",
    "text": ""
},

{
    "location": "lib/transformations.html#Transformations-1",
    "page": "Transformations",
    "title": "Transformations",
    "category": "section",
    "text": "This module provides functions to apply coordinate transformations to Systems using matrix decompositions.Pages = [\"transformations.md\"]\nDepth = 3CurrentModule = Reachability.Transformations"
},

{
    "location": "lib/transformations.html#Reachability.Transformations.transform",
    "page": "Transformations",
    "title": "Reachability.Transformations.transform",
    "category": "function",
    "text": "transform(S; [method])\n\nInterface function that calls the respective transformation function.\n\nInput\n\nS      – discrete or continuous system\nmethod – (optional, default: \'schur\') transformation method name; valid             otions are:\n\'schur\'\n\nOutput\n\nA tuple containing:\n\ntransformed discrete or continuous system\ninverse transformation matrix for reverting the transformation\n\nNotes\n\nThe functions that are called in the background should return a the transformed system components A, X0, and U, and also an inverse transformation matrix M.\n\n\n\n"
},

{
    "location": "lib/transformations.html#Interface-1",
    "page": "Transformations",
    "title": "Interface",
    "category": "section",
    "text": "This module exports a single function that works as an interface. It dispatches which transformation to apply using a string argument.transform"
},

{
    "location": "lib/transformations.html#Reachability.Transformations.schur_transform",
    "page": "Transformations",
    "title": "Reachability.Transformations.schur_transform",
    "category": "function",
    "text": "schur_transform(S)\n\nApplies a Schur transformation to a discrete or continuous system.\n\nInput\n\nS – discrete or continuous system\n\nOutput\n\nA tuple containing:\n\ntransformed discrete or continuous system\ninverse transformation matrix for reverting the transformation\n\nAlgorithm\n\nWe use Julia\'s default schurfact function to compute a Schur decomposition of the coefficients matrix A.\n\n\n\n"
},

{
    "location": "lib/transformations.html#Schur-transform-1",
    "page": "Transformations",
    "title": "Schur transform",
    "category": "section",
    "text": "The real Schur decomposition is of the formU^TAU = beginpmatrix\nT_11  T_12 cdots  T_1b \n0  T_22  cdots  T_2b \nvdots  vdots  ddots  vdots \n0  0  cdots  T_bb\nendpmatrixwhere T_ij are 2x2 matricesschur_transform"
},

{
    "location": "lib/transformations.html#Examples-1",
    "page": "Transformations",
    "title": "Examples",
    "category": "section",
    "text": ""
},

{
    "location": "lib/discretize.html#",
    "page": "Discretization",
    "title": "Discretization",
    "category": "page",
    "text": ""
},

{
    "location": "lib/discretize.html#Reachability.ReachSets.discretize",
    "page": "Discretization",
    "title": "Reachability.ReachSets.discretize",
    "category": "function",
    "text": "discretize(cont_sys, δ; [approx_model], [pade_expm], [lazy_expm], [lazy_sih])\n\nDiscretize a continuous system of ODEs with nondeterministic inputs.\n\nInput\n\ncont_sys     – continuous system\nδ            – step size\napprox_model – the method to compute the approximation model for the                   discretization, among:\nforward    – use forward-time interpolation\nbackward   – use backward-time interpolation\nfirstorder – use first order approximation of the ODE\nnobloating – do not bloat the initial states                 (use for discrete-time reachability)\npade_expm    – (optional, default = false) if true, use Pade approximant                   method to compute matrix exponentials of sparse matrices;                   otherwise use Julia\'s buil-in expm\nlazy_expm    – (optional, default = false) if true, compute the matrix                   exponential in a lazy way (suitable for very large systems)\nlazy_sih     – (optional, default = true) if true, compute the                   symmetric interval hull in a lazy way (suitable if only a                   few dimensions are of interest)\n\nOutput\n\nA discrete system.\n\nNotes\n\nThis function applies an approximation model to transform a continuous affine system into a discrete affine system. This transformation allows to do dense time reachability, i.e. such that the trajectories of the given continuous system are included in the computed flowpipe of the discretized system. For discrete-time reachability, use approx_model=\"nobloating\".\n\n\n\n"
},

{
    "location": "lib/discretize.html#Discretize-1",
    "page": "Discretization",
    "title": "Discretize",
    "category": "section",
    "text": "Pages = [\"discretize.md\"]\nDepth = 3CurrentModule = Reachability.ReachSetsdiscretize"
},

{
    "location": "lib/discretize.html#Reachability.ReachSets.discr_bloat_firstorder",
    "page": "Discretization",
    "title": "Reachability.ReachSets.discr_bloat_firstorder",
    "category": "function",
    "text": "bloat_firstorder(cont_sys, δ)\n\nCompute bloating factors using first order approximation.\n\nInput\n\ncont_sys – a continuous affine system\nδ        – step size\n\nNotes\n\nIn this algorithm, the infinity norm is used. See also: discr_bloat_interpolation for more accurate (less conservative) bounds.\n\nAlgorithm\n\nThis uses a first order approximation of the ODE, and matrix norm upper bounds, see Le Guernic, C., & Girard, A., 2010, Reachability analysis of linear systems using support functions. Nonlinear Analysis: Hybrid Systems, 4(2), 250-262.\n\n\n\n"
},

{
    "location": "lib/discretize.html#Reachability.ReachSets.discr_bloat_interpolation",
    "page": "Discretization",
    "title": "Reachability.ReachSets.discr_bloat_interpolation",
    "category": "function",
    "text": "discr_bloat_interpolation(cont_sys, δ, approx_model, pade_expm, lazy_expm)\n\nCompute bloating factors using forward or backward interpolation.\n\nInput\n\ncs           – a continuous system\nδ            – step size\napprox_model – choose the approximation model among \"forward\" and                   \"backward\"\npade_expm    – if true, use Pade approximant method to compute the                   matrix exponential\nlazy_expm   –  if true, compute the matrix exponential in a lazy way                   suitable for very large systems)\n\nAlgorithm\n\nSee Frehse et al., CAV\'11, SpaceEx: Scalable Verification of Hybrid Systems, Lemma 3.\n\nNote that in the unlikely case that A is invertible, the result can also be obtained directly, as a function of the inverse of A and e^{At} - I.\n\nThe matrix P is such that: ϕAabs = P[1:n, 1:n], Phi1Aabsdelta = P[1:n, (n+1):2*n], and Phi2Aabs = P[1:n, (2*n+1):3*n].\n\n\n\n"
},

{
    "location": "lib/discretize.html#Dense-time-reachability-1",
    "page": "Discretization",
    "title": "Dense-time reachability",
    "category": "section",
    "text": "discr_bloat_firstorder\ndiscr_bloat_interpolation"
},

{
    "location": "lib/discretize.html#Reachability.ReachSets.discr_no_bloat",
    "page": "Discretization",
    "title": "Reachability.ReachSets.discr_no_bloat",
    "category": "function",
    "text": "discr_no_bloat(cont_sys, δ, pade_expm, lazy_expm)\n\nDiscretize a continuous system without bloating of the initial states, suitable for discrete-time reachability.\n\nInput\n\ncont_sys     – a continuous system\nδ            – step size\npade_expm    – if true, use Pade approximant method to compute the                   matrix exponential\nlazy_expm    – if true, compute the matrix exponential in a lazy way                   (suitable for very large systems)\n\nOutput\n\nA discrete system.\n\nAlgorithm\n\nThe transformation implemented here is the following:\n\nA -> Phi := exp(A*delta)\nU -> V := M*U\nX0 -> X0hat := X0\n\nwhere M corresponds to Phi1(A, delta) in Eq. (8) of SpaceEx: Scalable Verification of Hybrid Systems.\n\nIn particular, there is no bloating, i.e. we don\'t bloat the initial states and dont multiply the input by the step size δ, as required for the dense time case.\n\n\n\n"
},

{
    "location": "lib/discretize.html#Discrete-time-reachability-1",
    "page": "Discretization",
    "title": "Discrete-time reachability",
    "category": "section",
    "text": "discr_no_bloat"
},

{
    "location": "lib/distributed.html#",
    "page": "Distributed computations",
    "title": "Distributed computations",
    "category": "page",
    "text": ""
},

{
    "location": "lib/distributed.html#Distributed-computations-1",
    "page": "Distributed computations",
    "title": "Distributed computations",
    "category": "section",
    "text": "This section of the manual describes functions to make use of distributed computation.Pages = [\"distributed.md\"]CurrentModule = Reachability"
},

{
    "location": "lib/distributed.html#Using-multiple-threads-1",
    "page": "Distributed computations",
    "title": "Using multiple threads",
    "category": "section",
    "text": "To control the number of threads used by your BLAS library, use the function Base.LinAlg.BLAS.set_num_threads(n), where n is an integer. Furthermore, the function get_num_threads() defined below will return the current value.Note. If you are using Julia v\"0.7-\" (run the command VERSION to find this), instead of Base.LinAlg below use LinearAlgebra, and this module should have been loaded in the current scope with using LinearAlgebra.#\n# This function is a part of Julia. License is MIT: https://julialang.org/license\n#\nfunction get_num_threads() # anonymous so it will be serialized when called\n    blas = Base.LinAlg.BLAS.vendor()\n    # Wrap in a try to catch unsupported blas versions\n    try\n        if blas == :openblas\n            return ccall((:openblas_get_num_threads, Base.libblas_name), Cint, ())\n        elseif blas == :openblas64\n            return ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())\n        elseif blas == :mkl\n            return ccall((:MKL_Get_Max_Num_Threads, Base.libblas_name), Cint, ())\n        end\n\n        # OSX BLAS looks at an environment variable\n        if Sys.isapple()\n            return ENV[\"VECLIB_MAXIMUM_THREADS\"]\n        end\n    end\n\n    return nothing\nend"
},

{
    "location": "about.html#",
    "page": "About",
    "title": "About",
    "category": "page",
    "text": ""
},

{
    "location": "about.html#About-1",
    "page": "About",
    "title": "About",
    "category": "section",
    "text": "This page contains some general information about this project, and recommendations about contributing.Pages = [\"about.md\"]"
},

{
    "location": "about.html#Contributing-1",
    "page": "About",
    "title": "Contributing",
    "category": "section",
    "text": "If you like this package, consider contributing! We welcome bug reports, examples that can be added to the documentation, or new feature proposals.Below some conventions that we follow when contributing to this package are detailed. For specific guidelines on documentation, see the Documentations Guidelines wiki."
},

{
    "location": "about.html#Branches-1",
    "page": "About",
    "title": "Branches",
    "category": "section",
    "text": "Each pull request (PR) should be pushed in a new branch with the name of the author followed by a descriptive name, e.g. t/mforets/my_feature. If the branch is associated to a previous discussion in one issue, we use the name of the issue for easier lookup, e.g. t/mforets/7."
},

{
    "location": "about.html#Unit-testing-and-continuous-integration-(CI)-1",
    "page": "About",
    "title": "Unit testing and continuous integration (CI)",
    "category": "section",
    "text": "This project is synchronized with Travis CI, such that each PR gets tested before merging (and the build is automatically triggered after each new commit). For the maintainability of this project, it is important to understand and fix the failing doctests if they exist. We develop in Julia v0.6.0, but for experimentation we also build on the nightly branch.To run the unit tests locally, you should do:$ julia --color=yes test/runtests.jl"
},

{
    "location": "about.html#Contributing-to-the-documentation-1",
    "page": "About",
    "title": "Contributing to the documentation",
    "category": "section",
    "text": "This documentation is written in Markdown, and it relies on Documenter.jl to produce the HTML layout. To build the docs, run make.jl:$ julia --color=yes docs/make.jl"
},

{
    "location": "about.html#References-1",
    "page": "About",
    "title": "References",
    "category": "section",
    "text": "This repository was originally motivated by the mathematical approach described in Reach Set Approximation through Decomposition with Low-dimensional Sets and High-dimensional Matrices,  Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Andreas Podelski, Christian Schilling, Frédéric Viry, in 21st ACM International Conference on Hybrid Systems: Computation and Control, 2018 Edition (Porto, Portugal), see the arXiv pre-print here.For a full references list of algorithms implemented in this repository, consult the References wiki."
},

{
    "location": "about.html#Credits-1",
    "page": "About",
    "title": "Credits",
    "category": "section",
    "text": "These persons have contributed to Reachability.jl (in alphabetic order):Marcelo Forets\nChristian Schilling\nFrederic ViryWe are also grateful to Goran Frehse, Sergiy Bogomolov, Andreas Podelski, Alexandre Rocca and Nikolaos Kekatos for enlightening discussions."
},

]}
