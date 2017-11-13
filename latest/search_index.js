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
    "text": "Pages = [\n    \"lib/interface.md\",\n    \"lib/discretize.md\"\n]\nDepth = 2"
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
    "category": "Function",
    "text": "solve(system, options)  or  solve(system, :key1 => val1, [...], keyK => valK)\n\nSolves a reachability problem s.t. the given options. If some options are not defined, we may fall back to default values.\n\nInput\n\nsystem  – a (discrete or continuoues) system specification\noptions – options for solving the problem\n\n\n\n"
},

{
    "location": "lib/interface.html#Posing-and-solving-a-reachability-problem-1",
    "page": "User interface",
    "title": "Posing and solving a reachability problem",
    "category": "section",
    "text": "solve"
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
    "category": "Function",
    "text": "discretize(cont_sys, δ; [approx_model], [pade_expm], [lazy_expm])\n\nDiscretize a continuous system of ODEs with nondeterministic inputs.\n\nInput\n\ncont_sys          – continuous system\nδ                 – step size\napprox_model      – the method to compute the approximation model for the                        discretization, among:\nforward    – use forward-time interpolation\nbackward   – use backward-time interpolation\nfirstorder – use first order approximation of the ODE\nnobloating – do not bloat the initial states                 (use for discrete-time reachability)\npade_expm         – (optional, default = false) if true, use Pade approximant                        method to compute matrix exponentials of sparse matrices;                        otherwise use Julia's buil-in expm\nlazy_expm         – (optional, default = false) if true, compute the matrix                        exponential in a lazy way (suitable for very large systems)\n\nOutput\n\nA discrete system.\n\nNotes\n\nThis function applies an approximation model to transform a continuous affine system into a discrete affine system. This transformation allows to do dense time reachability, i.e. such that the trajectories of the given continuous system are included in the computed flowpipe of the discretized system. For discrete-time reachability, use approx_model=\"nobloating\".\n\n\n\n"
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
    "category": "Function",
    "text": "bloat_firstorder(cont_sys, δ)\n\nCompute bloating factors using first order approximation.\n\nInput\n\ncont_sys – a continuous affine system\nδ        – step size\n\nNotes\n\nIn this algorithm, the infinity norm is used. See also: discr_bloat_interpolation for more accurate (less conservative) bounds.\n\nAlgorithm\n\nThis uses a first order approximation of the ODE, and matrix norm upper bounds, see Le Guernic, C., & Girard, A., 2010, Reachability analysis of linear systems using support functions. Nonlinear Analysis: Hybrid Systems, 4(2), 250-262.\n\n\n\n"
},

{
    "location": "lib/discretize.html#Reachability.ReachSets.discr_bloat_interpolation",
    "page": "Discretization",
    "title": "Reachability.ReachSets.discr_bloat_interpolation",
    "category": "Function",
    "text": "discr_bloat_interpolation(cont_sys, δ, approx_model, pade_expm, lazy_expm)\n\nCompute bloating factors using forward or backward interpolation.\n\nInput\n\ncs        – a continuous system\nδ         – step size\napprox_model\npade_expm – if true, use Pade approximant method to compute the                matrix exponential\nlazy_expm – if true, compute the matrix exponential in a lazy way                suitable for very large systems)\n\nAlgorithm\n\nSee Frehse et al CAV'11 paper, SpaceEx: Scalable Verification of Hybrid Systems, see Lemma 3.\n\nNote that in the unlikely case that A is invertible, the result can also be obtained directly, as a function of the inverse of A and e^{At} - I.\n\nThe matrix P is such that: ϕAabs = P[1:n, 1:n], Phi1Aabsdelta = P[1:n, (n+1):2*n], and Phi2Aabs = P[1:n, (2*n+1):3*n].\n\n\n\n"
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
    "category": "Function",
    "text": "discr_no_bloat(cont_sys, δ, pade_expm, lazy_expm)\n\nDiscretize a continuous system without bloating of the initial states, suitable for discrete-time reachability.\n\nInput\n\ncont_sys     – a continuous system\nδ            – step size\npade_expm    – if true, use Pade approximant method to compute the                   matrix exponential\nlazy_expm    – if true, compute the matrix exponential in a lazy way                   (suitable for very large systems)\n\nOutput\n\nA discrete system.\n\nAlgorithm\n\nThe transformation implemented here is the following:\n\nA -> Phi := exp(A*delta)\nU -> V := M*U\nX0 -> X0hat := X0\n\nwhere M corresponds to Phi1(A, delta) in Eq. (8) of SpaceEx: Scalable Verification of Hybrid Systems.\n\nIn particular, there is no bloating, i.e. we don't bloat the initial states and dont multiply the input by the step size δ, as required for the dense time case.\n\n\n\n"
},

{
    "location": "lib/discretize.html#Discrete-time-reachability-1",
    "page": "Discretization",
    "title": "Discrete-time reachability",
    "category": "section",
    "text": "discr_no_bloat"
},

{
    "location": "about.html#",
    "page": "About",
    "title": "About",
    "category": "page",
    "text": ""
},

{
    "location": "about.html#Contributing-1",
    "page": "About",
    "title": "Contributing",
    "category": "section",
    "text": "Pages = [\"about/CONTRIBUTING.md\"]This page details the some of the guidelines that should be followed when contributing to this package."
},

{
    "location": "about.html#Running-the-Unit-Tests-1",
    "page": "About",
    "title": "Running the Unit Tests",
    "category": "section",
    "text": "$ julia --color=yes test/runtests.jl"
},

{
    "location": "about.html#Branches-1",
    "page": "About",
    "title": "Branches",
    "category": "section",
    "text": ""
},

{
    "location": "about.html#Contributing-to-the-Documentation-1",
    "page": "About",
    "title": "Contributing to the Documentation",
    "category": "section",
    "text": "The documentation source is written with Markdown, and we use Documenter.jl to produce the HTML documentation. To build the docs, run make.jl:$ julia --color=yes docs/make.jl"
},

{
    "location": "about.html#Credits-1",
    "page": "About",
    "title": "Credits",
    "category": "section",
    "text": "These persons have contributed to Reachability.jl (in alphabetic order):Marcelo Forets\nChristian Schilling\nFrederic ViryWe are also grateful to Goran Frehse for enlightening discussions."
},

]}
