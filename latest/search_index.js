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
    "text": "Pages = [\n    \"lib/interface.md\",\n]\nDepth = 2"
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
    "text": "solve(system, options)  or  solve(system, :key1 => val1, [...], keyK => valK)\n\nSolves a reachability problem s.t. the given options. If some options are not defined, we may fall back to default values.\n\nINPUT:\n\nsystem – a (discrete or continuoues) system specification\noptions – options for solving the problem\n\n\n\n"
},

{
    "location": "lib/interface.html#Posing-and-solving-a-reachability-problem-1",
    "page": "User interface",
    "title": "Posing and solving a reachability problem",
    "category": "section",
    "text": "solve"
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
