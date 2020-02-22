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
    "text": "Pages = [\n    \"lib/interface.md\",\n    \"lib/systems.md\",\n    \"lib/algorithms.md\",\n    \"lib/transformations.md\",\n    \"lib/discretize.md\",\n    \"lib/distributed.md\"\n]\nDepth = 2"
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
    "text": "solve(system, options)  or  solve(system, :key1 => val1, [...], keyK => valK)\n\nSolves a reachability problem s.t. the given options. If some options are not defined, we may fall back to default values.\n\nInput\n\nsystem    – a (discrete or continuoues) system specification\noptions   – algorithm options for solving the problem\nalgorithm – (optional, default: dispatched on the system\'s type) the                reachability algorithm for the computation\n\nOutput\n\nA solution object whose content depends on the input options.\n\nNotes\n\nTo see all available input options, see keys(Reachability.available_keywords.dict).\n\n\n\n\n\nsolve(system::InitialValueProblem{<:HybridSystem},\n      options::Options)::AbstractSolution\n\nInterface to reachability algorithms for a hybrid system PWA dynamics.\n\nInput\n\nsystem  – hybrid system\noptions – options for solving the problem\n\n\n\n\n\n"
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
    "text": "This module provides types to represent systems of affine ODEs in both discrete and continuous time.Pages = [\"systems.md\"]\nDepth = 3"
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
    "location": "lib/systems.html#Reachability.Utils.normalize",
    "page": "Systems",
    "title": "Reachability.Utils.normalize",
    "category": "function",
    "text": "normalize(system::AbstractSystem)\n\nTransform a mathematical system to a normalized (or canonical) form.\n\nInput\n\nsystem – system; it can be discrete or continuous\n\nOutput\n\nEither the same system if it already conforms to a canonical form, or a new system otherwise.\n\nNotes\n\nThe normalization procedure consists of transforming a given system type into a \"canonical\" format that is used internally. More details are given below.\n\nAlgorithm\n\nThe implementation of normalize exploits MathematicalSystems\'s\' types, which carry information about the problem as a type parameter.\n\nHomogeneous ODEs of the form x = Ax x  mathcalX are canonical if the associated problem is a ConstrainedLinearContinuousSystem and A is a matrix. This type does not handle non-deterministic inputs.\n\nNote that a LinearContinuousSystem does not consider constraints on the state-space (such as an invariant); to specify state constraints, use a ConstrainedLinearContinuousSystem. If the passed system is a LinearContinuousSystem (i.e. no constraints) then the normalization fixes a universal set (Universe) as the constraint set.\n\nThe generalization to canonical systems with constraints and possibly time-varying non-deterministic inputs is considered next. These systems are of the form x = Ax + u u  mathcalU x  mathcalX. The system type is ConstrainedLinearControlContinuousSystem, where A is a matrix, X is a set and U is an input, that is, any concrete subtype of AbstractInput.\n\nIf U is not given as an input, normalization accepts either a LazySet, or a vector of LazySets. In these cases, the sets are wrapped around an appropriate concrete input type.\n\nIf the system does not conform to a canonical form, the implementation tries to make the transformation; otherwise an error is thrown. In particular, ODEs of the form x = Ax + Bu are mapped into x = Ax + u u  BmathcalU, where now u has the same dimensions as x.\n\nThe transformations described above are analogous in the discrete case, i.e. x_k+1 = A x_k and x_k+1 = Ax_k + u_k u_k  mathcalU x_k  mathcalX for the linear and affine cases respectively.\n\n\n\n\n\n"
},

{
    "location": "lib/systems.html#Normalization-1",
    "page": "Systems",
    "title": "Normalization",
    "category": "section",
    "text": "normalize"
},

{
    "location": "lib/algorithms.html#",
    "page": "Algorithms",
    "title": "Algorithms",
    "category": "page",
    "text": ""
},

{
    "location": "lib/algorithms.html#Available-Algorithms-1",
    "page": "Algorithms",
    "title": "Available Algorithms",
    "category": "section",
    "text": "This section of the manual describes the algorithms that are available in this package.Pages = [\"interface.md\"]CurrentModule = Reachability.ReachSets"
},

{
    "location": "lib/algorithms.html#Continuous-time-reachability-1",
    "page": "Algorithms",
    "title": "Continuous-time reachability",
    "category": "section",
    "text": ""
},

{
    "location": "lib/algorithms.html#Reachability.ReachSets.BFFPSV18",
    "page": "Algorithms",
    "title": "Reachability.ReachSets.BFFPSV18",
    "category": "type",
    "text": "BFFPSV18 <: ContinuousPost\n\nImplementation of the reachability algorithm for purely continuous linear time-invariant systems using block decompositons by S. Bogomolov, M. Forets, G. Frehse, A. Podelski, C. Schilling and F. Viry [1].\n\nFields\n\noptions – an Options structure that holds the algorithm-specific options\n\nNotes\n\nThe following options are available:\n\noption :discretization of type String has default value \'forward\'; model for bloating/continuous time analysis\noption :algorithm of type String has default value \'explicit\'; algorithm backend\noption :δ of type Float64 with alias :sampling_time has default value \'0.01\'; time step\noption :vars of type AbstractArray{Int64,1} has default value \'Int64[]\'; variables of interest; default: all variables\noption :partition of type AbstractArray{#s124,1} where #s124<:AbstractArray{Int64,1} has default value \'Array{Int64,1}[[]]\'; block partition; a block is represented by a vector containing its indices\noption :sih_method of type String has default value \'concrete\'; method to compute the symmetric interval hull in discretization\noption :exp_method of type String has default value \'base\'; method to compute the matrix exponential\noption :assume_sparse of type Bool has default value \'false\'; use an analysis for sparse discretized matrices?\noption :lazy_inputs_interval of type Union{Int64, Function} has default value \'getfield(Reachability.ReachSets, Symbol(\"##24#25\"))()\'; length of interval in which the inputs are handled as a lazy set (``-1`` for \'never\'); may generally also be a predicate over indices; the default corresponds to ``-1``\noption :block_options of type Any has default value \'nothing\'; short hand to set \':block_options_init\' and \':block_options_iter\'\noption :block_options_init of type Any has default value \'nothing\'; option for the approximation of the initial states (during decomposition)\noption :block_options_iter of type Any has default value \'nothing\'; option for the approximation of the states ``X_k``, ``k>0``\noption :assume_homogeneous of type Bool has default value \'false\'; ignore dynamic inputs during the analysis?\noption :eager_checking of type Bool has default value \'true\'; terminate as soon as property violation was detected?\n\n\nAlgorithm\n\nWe refer to [1] for technical details.\n\n[1] Reach Set Approximation through Decomposition with Low-dimensional Sets and High-dimensional Matrices. S. Bogomolov, M. Forets, G. Frehse, A. Podelski, C. Schilling, F. Viry. HSCC \'18 Proceedings of the 21st International Conference on Hybrid Systems: Computation and Control (part of CPS Week).\n\n\n\n\n\n"
},

{
    "location": "lib/algorithms.html#Decomposition-based-approach-1",
    "page": "Algorithms",
    "title": "Decomposition-based approach",
    "category": "section",
    "text": "BFFPSV18"
},

{
    "location": "lib/properties.html#",
    "page": "Properties",
    "title": "Properties",
    "category": "page",
    "text": ""
},

{
    "location": "lib/properties.html#Properties-1",
    "page": "Properties",
    "title": "Properties",
    "category": "section",
    "text": "This module provides representations of (safety) properties.Pages = [\"properties.md\"]\nDepth = 3CurrentModule = Reachability.Properties"
},

{
    "location": "lib/properties.html#Reachability.Properties.Property",
    "page": "Properties",
    "title": "Reachability.Properties.Property",
    "category": "type",
    "text": "Abstract supertype of properties that can be checked.\n\nEvery concrete subtype should provide the following functions:\n\ndim(𝑃::Property)::Int\ncheck(𝑃::Property, X::LazySet; witness::Bool=false)\n\n\n\n\n\n"
},

{
    "location": "lib/properties.html#General-property-interface-1",
    "page": "Properties",
    "title": "General property interface",
    "category": "section",
    "text": "Property"
},

{
    "location": "lib/properties.html#Reachability.Properties.Conjunction",
    "page": "Properties",
    "title": "Reachability.Properties.Conjunction",
    "category": "type",
    "text": "Conjunction <: Property\n\nType that represents a conjunction of properties.\n\nFields\n\nconjuncts – vector of properties\n\nNotes\n\nThe following formula characterizes whether a set X satisfies a disjunction 𝑃 = 𝑃_1  𝑃_2    𝑃_m:\n\n    X models 𝑃 iff X models 𝑃_j text for all  1  j  m\n\n\n\n\n\n"
},

{
    "location": "lib/properties.html#LazySets.dim-Tuple{Conjunction}",
    "page": "Properties",
    "title": "LazySets.dim",
    "category": "method",
    "text": "dim(𝑃::Conjunction)::Int\n\nReturn the dimension of a conjunction of properties.\n\nInput\n\n𝑃 – conjunction of properties\n\nOutput\n\nThe dimension of the conjunction of properties.\n\n\n\n\n\n"
},

{
    "location": "lib/properties.html#Reachability.Properties.check-Union{Tuple{N}, Tuple{Conjunction,LazySet{N}}} where N<:Real",
    "page": "Properties",
    "title": "Reachability.Properties.check",
    "category": "method",
    "text": "check(𝑃::Conjunction, X::LazySet{N}; witness::Bool=false) where {N<:Real}\n\nCheck whether a convex set satisfies a conjunction of properties.\n\nInput\n\n𝑃       – conjunction of properties\nX       – convex set\nwitness – (optional, default: false) flag for returning a counterexample              if the property is violated\n\nOutput\n\nIf witness option is deactivated: true iff X satisfies the property 𝑃\nIf witness option is activated:\n(true, []) iff X satisfies the property 𝑃\n(false, v) iff X does not satisfy the property 𝑃 with witness v\n\nNotes\n\nBy convention, the empty conjunction is equivalent to true and hence is satisfied by any set.\n\n\n\n\n\n"
},

{
    "location": "lib/properties.html#Reachability.Properties.Disjunction",
    "page": "Properties",
    "title": "Reachability.Properties.Disjunction",
    "category": "type",
    "text": "Disjunction <: Property\n\nType that represents a disjunction of properties.\n\nFields\n\ndisjuncts – vector of properties (elements are reordered by this type)\nreorder   – flag to indicate whether shuffling is allowed\n\nNotes\n\nThe following formula characterizes whether a set X satisfies a disjunction 𝑃 = 𝑃_1  𝑃_2    𝑃_m:\n\n    X models 𝑃 iff X models 𝑃_j text for some  1  j  m\n\nIf the reorder flag is set, the disjuncts may be reordered after each call to check as a heuristics to make subsequent checks faster.\n\n\n\n\n\n"
},

{
    "location": "lib/properties.html#LazySets.dim-Tuple{Disjunction}",
    "page": "Properties",
    "title": "LazySets.dim",
    "category": "method",
    "text": "dim(𝑃::Disjunction)::Int\n\nReturn the dimension of a disjunction of properties.\n\nInput\n\n𝑃 – disjunction of properties\n\nOutput\n\nThe dimension of the disjunction of properties.\n\n\n\n\n\n"
},

{
    "location": "lib/properties.html#Reachability.Properties.check-Union{Tuple{N}, Tuple{Disjunction,LazySet{N}}} where N<:Real",
    "page": "Properties",
    "title": "Reachability.Properties.check",
    "category": "method",
    "text": "check(𝑃::Disjunction, X::LazySet{N}; witness::Bool=false) where {N<:Real}\n\nCheck whether a convex set satisfies a disjunction of properties.\n\nInput\n\n𝑃       – disjunction of properties\nX       – convex set\nwitness – (optional, default: false) flag for returning a counterexample              if the property is violated\n\nOutput\n\nIf witness option is deactivated: true iff X satisfies the property 𝑃\nIf witness option is activated:\n(true, []) iff X satisfies the property 𝑃\n(false, v) iff X does not satisfy the property 𝑃 with witness v; note that v == N[] if 𝑃 is the empty disjunction\n\nNotes\n\nBy convention, the empty disjunction is equivalent to false and hence is satisfied by no set.\n\nIf the 𝑃.reorder flag is set, the disjuncts may be reordered as a heuristics to make subsequent checks faster. Since we check satisfaction from left to right, we move the disjunct for which satisfaction was established to the front.\n\nTo be consistent with other propertes, the witness option only returns one counterexample, namely for the left-most disjunct in the disjuncts vector. We deactivate witness production for checking the remaining disjuncts.\n\n\n\n\n\n"
},

{
    "location": "lib/properties.html#Boolean-combination-of-properties-1",
    "page": "Properties",
    "title": "Boolean combination of properties",
    "category": "section",
    "text": "Conjunction\ndim(::Conjunction)\ncheck(::Conjunction, ::LazySet{N}) where {N<:Real}\nDisjunction\ndim(::Disjunction)\ncheck(::Disjunction, ::LazySet{N}) where {N<:Real}"
},

{
    "location": "lib/properties.html#Reachability.Properties.SafeStatesProperty",
    "page": "Properties",
    "title": "Reachability.Properties.SafeStatesProperty",
    "category": "type",
    "text": "SafeStatesProperty{N<:Real} <: Property\n\nType that represents a safety property characterized by a set of safe states. The property is satisfied by a given set of states X if X is fully contained in the set of safe states.\n\nFields\n\nsafe    – convex set representing the safe states\nwitness – witness point (empty vector if not set)\n\nNotes\n\nThe following formula characterizes whether a set X satisfies a safety property characterized by a set of safe states 𝑃:\n\n    X models 𝑃 iff X  𝑃textttsafe\n\n\n\n\n\n"
},

{
    "location": "lib/properties.html#LazySets.dim-Tuple{SafeStatesProperty}",
    "page": "Properties",
    "title": "LazySets.dim",
    "category": "method",
    "text": "dim(𝑃::SafeStatesProperty)::Int\n\nReturn the dimension of a property with safe states.\n\nInput\n\n𝑃 – safety property with safe states\n\nOutput\n\nThe dimension of the safe states.\n\n\n\n\n\n"
},

{
    "location": "lib/properties.html#Reachability.Properties.check-Tuple{SafeStatesProperty,LazySet}",
    "page": "Properties",
    "title": "Reachability.Properties.check",
    "category": "method",
    "text": "check(𝑃::SafeStatesProperty, X::LazySet; witness::Bool=false)\n\nChecks whether a convex set is contained in the set of safe states.\n\nInput\n\n𝑃       – safety property with safe states\nX       – convex set\nwitness – (optional, default: false) flag for returning a counterexample              if the property is violated\n\nOutput\n\nLet Y be the safe states represented by 𝑃.\n\nIf witness option is deactivated: true iff X  Y\nIf witness option is activated:\n(true, []) iff X  Y\n(false, v) iff X  Y and v  X setminus Y\n\n\n\n\n\n"
},

{
    "location": "lib/properties.html#Reachability.Properties.BadStatesProperty",
    "page": "Properties",
    "title": "Reachability.Properties.BadStatesProperty",
    "category": "type",
    "text": "BadStatesProperty{N<:Real} <: Property\n\nType that represents a safety property characterized by a set of bad states. The property is satisfied by a given set of states if the intersection with the set of bad states is empty.\n\nFields\n\nbad     – convex set representing the bad states\nwitness – witness point (empty vector if not set)\n\nNotes\n\nThe following formula characterizes whether a set X satisfies a safety property characterized by a set of bad states 𝑃:\n\n    X models 𝑃 iff X  𝑃textttbad = \n\n\n\n\n\n"
},

{
    "location": "lib/properties.html#LazySets.dim-Tuple{BadStatesProperty}",
    "page": "Properties",
    "title": "LazySets.dim",
    "category": "method",
    "text": "dim(𝑃::BadStatesProperty)::Int\n\nReturn the dimension of a property with bad states.\n\nInput\n\n𝑃 – safety property with bad states\n\nOutput\n\nThe dimension of the bad states.\n\n\n\n\n\n"
},

{
    "location": "lib/properties.html#Reachability.Properties.check-Tuple{BadStatesProperty,LazySet}",
    "page": "Properties",
    "title": "Reachability.Properties.check",
    "category": "method",
    "text": "check(𝑃::BadStatesProperty, X::LazySet; witness::Bool=false)\n\nChecks whether a convex set is disjoint from the set of bad states.\n\nInput\n\n𝑃       – safety property with bad states\nX       – convex set\nwitness – (optional, default: false) flag for returning a counterexample              if the property is violated\n\nOutput\n\nLet Y be the bad states represented by 𝑃.\n\nIf witness option is deactivated: true iff X  Y = \nIf witness option is activated:\n(true, []) iff X  Y = \n(false, v) iff X  Y   and v  X  Y\n\n\n\n\n\n"
},

{
    "location": "lib/properties.html#Specific-properties-1",
    "page": "Properties",
    "title": "Specific properties",
    "category": "section",
    "text": "SafeStatesProperty\ndim(::SafeStatesProperty)\ncheck(::SafeStatesProperty, ::LazySet)\nBadStatesProperty\ndim(::BadStatesProperty)\ncheck(::BadStatesProperty, ::LazySet)"
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
    "location": "lib/transformations.html#Reachability.Utils.transform",
    "page": "Transformations",
    "title": "Reachability.Utils.transform",
    "category": "function",
    "text": "transform(problem::InitialValueProblem, options::Options)\n\nInterface function that calls the respective transformation function.\n\nInput\n\nproblem – discrete or continuous initial-value problem\noption  – problem options\n\nOutput\n\nA tuple containing the transformed problem and the transformed options.\n\nNotes\n\nThe functions that are called in the background should return a the transformed system components A, X0, and U, and also an inverse transformation matrix M. If the system has an invariant, it is transformed as well.\n\n\n\n\n\n"
},

{
    "location": "lib/transformations.html#Interface-1",
    "page": "Transformations",
    "title": "Interface",
    "category": "section",
    "text": "This module exports a single function that works as an interface. It dispatches which transformation to apply using a string argument.transform"
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
    "text": "discretize(𝑆, δ; [algorithm], [exp_method], [sih_method], [set_operations])\n\nApply an approximation model to S obtaining a discrete initial value problem.\n\nInput\n\n𝑆             – initial value problem for a continuous affine ODE with                    non-deterministic inputs\nδ             – step size\nalgorithm     – (optional, default: \"forward\") the algorithm used to                    compute the approximation model for the discretization,                    choose among:\n\"forward\"    – use forward-time interpolation\n\"backward\"   – use backward-time interpolation\n\"firstorder\" – use first-order approximation of the ODE\n\"nobloating\" – do not bloat the initial states\nexp_method  – (optional, default: \"base\") the method used to take the matrix                   exponential of the coefficient matrix, choose among:\n\"base\" – the scaling and squaring method implemented in Julia base,             see ?exp for details\n\"pade\" – use Pade approximant method to compute matrix exponentials of             sparse matrices, implemented in Expokit\n\"lazy\" – compute a wrapper type around the matrix exponential, i.e. using             the lazy implementation SparseMatrixExp from LazySets and             the evaluation of the action of the matrix exponential using the             expmv implementation from Expokit\nsih_method  – (optional, default: \"concrete\") the method used to take the                   symmetric interval hull operation, choose among:\n\"concrete\" – compute the full symmetric interval hull using the function                 symmetric_interval_hull from LazySets.Approximations\n\"lazy\"     – compute a wrapper set type around symmetric interval hull                 in a lazy way using SymmetricIntervalHull\nset_operations  – (optional, default: \"lazy\") set operations used for the                      discretized initial states and transformed inputs, choose among:\n\"lazy\"     – use lazy convex hull for the initial states and lazy linear                 map for the inputs\n\"zonotope\" – use concrete zonotope operations (linear map and Minkowski sum)                 and return zonotopes for both the initial states and the                 inputs of the discretized system\n\nOutput\n\nThe initial value problem of a discrete system.\n\nAlgorithm\n\nLet 𝑆  x = Ax(t) + u(t), x(0)  mathcalX_0, u(t)  U be the given continuous affine ODE 𝑆, where U is the set of non-deterministic inputs and mathcalX_0 is the set of initial states. Recall that the system 𝑆 is called homogeneous whenever U is the empty set.\n\nGiven a step size δ, this function computes a set, Ω₀, that guarantees to contain all the trajectories of 𝑆 starting at any x(0)  mathcalX_0 and for any input function that satisfies u(t)  U, for any t  0 δ.\n\nThe initial value problem returned by this function consists of the set Ω₀ together with the coefficient matrix ϕ = e^Aδ and a transformed set of inputs if U is non-empty.\n\nIn the literature, the method to obtain Ω₀ is called the approximation model and different alternatives have been proposed. See the argument algorithm for available options. For the reference to the original papers, see the docstring of each method discretize_....\n\nIn the dense-time case, the transformation is such that the trajectories of the given continuous system are included in the computed flowpipe of the discretized system.\n\nIn the discrete-time case, there is no bloating of the initial states and the input is assumed to remain constant between sampled times. Use the option algorithm=\"nobloating\" for this setting.\n\nSeveral methods to compute the matrix exponential are availabe. Use exp_method to select one. For very large systems, computing the full matrix exponential is expensive hence it is preferable to compute the action of the matrix exponential over vectors when needed, e^{δA} v for each v. Use the option exp_method=\"lazy\" for this purpose.\n\n\n\n\n\n"
},

{
    "location": "lib/discretize.html#Reachability.ReachSets.discretize_interpolation",
    "page": "Discretization",
    "title": "Reachability.ReachSets.discretize_interpolation",
    "category": "function",
    "text": "discretize_interpolation(𝑆, δ; [algorithm], [exp_method], [sih_method])\n\nCompute bloating factors using forward or backward interpolation.\n\nInput\n\n𝑆             – a continuous system\nδ             – step size\nalgorithm     – choose the algorithm to compute the approximation model                    among \"forward\" and \"backward\"\nexp_method    – (optional, default: \"base\") the method used to take the matrix                    exponential of the coefficient matrix, choose among:\n\"base\" – the scaling and squaring method implemented in Julia base,             see ?exp for details\n\"pade\" – use Pade approximant method to compute matrix exponentials of             sparse matrices, implemented in Expokit\n\"lazy\" – compute a wrapper type around the matrix exponential, i.e. using             the lazy implementation SparseMatrixExp from LazySets and             the evaluation of the action of the matrix exponential using the             expmv implementation from Expokit\nsih_method    – (optional, default: \"concrete\") the method used to take the                    symmetric interval hull operation, choose among:\n\"concrete\" – compute the full symmetric interval hull\n\"lazy\"     – compute a wrapper set type around symmetric interval hull in a                 lazy way\n\nOutput\n\nThe initial value problem for a discrete system. In particular:\n\nif the input system is homogeneous, a LinearDiscreteSystem is returned,\notherwise a ConstrainedLinearControlDiscreteSystem is returned.\n\nAlgorithm\n\nLet us define some notation. Let\n\n𝑆  x = Ax(t) + u(t)\n\nwith x(0)  mathcalX_0, u(t)  U be the given continuous affine ODE 𝑆, where U is the set of non-deterministic inputs and mathcalX_0 is the set of initial states.\n\nThe transformations are:\n\nΦ  exp^Aδ,\nΩ₀  ConvexHull(mathcalX_0 ΦmathcalX_0  δU(0)  Eψ(U(0) δ)  E^+(mathcalX_0 δ)),\nV  δU(k)  Eψ(U(k) δ).\n\nHere we allow U to be a sequence of time varying non-deterministic input sets.\n\nFor the definition of the sets Eψ and E^+ see [1]. The  \"backward\" method uses E^-.\n\n[1] Frehse, Goran, et al. \"SpaceEx: Scalable verification of hybrid systems.\" International Conference on Computer Aided Verification. Springer, Berlin, Heidelberg, 2011.\n\n\n\n\n\n"
},

{
    "location": "lib/discretize.html#Reachability.ReachSets.discretize_firstorder",
    "page": "Discretization",
    "title": "Reachability.ReachSets.discretize_firstorder",
    "category": "function",
    "text": "discretize_firstorder(𝑆, δ; [p], [exp_method])\n\nApply a first-order approximation model to S obtaining a discrete initial value problem.\n\nInput\n\n𝑆           – initial value problem for a continuous affine ODE with                  non-deterministic inputs\nδ           – step size\np           – (optional, default: Inf) parameter in the considered norm\nexp_method  – (optional, default: \"base\") the method used to take the matrix                   exponential of the coefficient matrix, choose among:\n\"base\" – the scaling and squaring method implemented in Julia base,             see ?exp for details\n\"pade\" – use Pade approximant method to compute matrix exponentials of             sparse matrices, implemented in Expokit\n\"lazy\" – compute a wrapper type around the matrix exponential, i.e. using             the lazy implementation SparseMatrixExp from LazySets and             the evaluation of the action of the matrix exponential using the             expmv implementation from Expokit\n\nOutput\n\nThe initial value problem for a discrete system. In particular:\n\nif the input system is homogeneous, a LinearDiscreteSystem is returned,\notherwise a ConstrainedLinearControlDiscreteSystem is returned.\n\nAlgorithm\n\nLet us define some notation. Let\n\n𝑆  x = Ax(t) + u(t)\n\nwith x(0)  mathcalX_0, u(t)  U be the given continuous affine ODE 𝑆, where U is the set of non-deterministic inputs and mathcalX_0 is the set of initial states.\n\nDefine R_mathcalX_0 = max_x  mathcalX_0 x, D_{\\mathcal{X}_0} = \\max_{x, y ∈ \\mathcal{X}_0} ‖x-y‖andR_{U} = \\max_{u ∈ U} ‖u‖. If only the support functions of\\mathcal{X}_0andUare known, these values might be hard to compute for some norms. SeeNotes` below for details.\n\nLet Ω₀ be the set defined as:\n\nΩ₀ = ConvexHull(mathcalX_0 e^δAmathcalX_0  δU  αB_p)\n\nwhere α = (e^δ A - 1 - δA)*(R_mathcalX_0 + R_U  A) and B_p denotes the unit ball for the considered p-norm.\n\nIt is proved in [Lemma 1, 1] that the set of states reachable by S in the time interval 0 δ, which we denote R_0δ(mathcalX_0) here, is included in Ω₀:\n\nR_0δ(mathcalX_0)  Ω₀\n\nMoreover, if d_H(A B) denotes the Hausdorff distance between the sets A and B in mathbbR^n, then\n\nd_H(Ω₀ R_0δ(mathcalX_0))  frac14(e^δ A - 1) D_mathcalX_0 + 2α\n\nHence, the approximation error can be made arbitrarily small by choosing δ small enough.\n\nHere we allow U to be a sequence of time varying non-deterministic inputs.\n\nNotes\n\nIn this implementation, the infinity norm is used by default. Other usual norms are p=1 and p=2. However, note that not all norms are supported; see the documentation of ?norm in LazySets for the supported norms.\n\nSee also discretize_interpolation for an alternative algorithm that uses less conservative bounds.\n\n[1] Le Guernic, C., & Girard, A., 2010, Reachability analysis of linear systems using support functions. Nonlinear Analysis: Hybrid Systems, 4(2), 250-262.\n\n\n\n\n\n"
},

{
    "location": "lib/discretize.html#Reachability.ReachSets.discretize_nobloating",
    "page": "Discretization",
    "title": "Reachability.ReachSets.discretize_nobloating",
    "category": "function",
    "text": "discretize_nobloating(𝑆, δ; [exp_method])\n\nDiscretize a continuous system without bloating of the initial states, suitable for discrete-time reachability.\n\nInput\n\n𝑆          – a continuous system\nδ          – step size\nexp_method – (optional, default: \"base\") the method used to take the matrix                 exponential of the coefficient matrix, choose among:\n\"base\" – the scaling and squaring method implemented in Julia base,             see ?exp for details\n\"pade\" – use Pade approximant method to compute matrix exponentials of             sparse matrices, implemented in Expokit\n\"lazy\" – compute a wrapper type around the matrix exponential, i.e. using             the lazy implementation SparseMatrixExp from LazySets and             the evaluation of the action of the matrix exponential using the             expmv implementation from Expokit\n\nOutput\n\nThe initial value problem for a discrete system. In particular:\n\nif the input system is homogeneous, a LinearDiscreteSystem is returned,\notherwise a ConstrainedLinearControlDiscreteSystem is returned.\n\nAlgorithm\n\nLet us define some notation. Let\n\n𝑆  x = Ax(t) + u(t)\n\nwith x(0)  mathcalX_0, u(t)  U be the given continuous affine ODE 𝑆, where U is the set of non-deterministic inputs and mathcalX_0 is the set of initial states.\n\nThe approximation model implemented in this function is such that there is no bloating, i.e. we don\'t bloat the initial states and don\'t multiply the input by the step size δ, as required for the dense time case.\n\nThe transformations are:\n\nΦ  exp^Aδ\nΩ₀  mathcalX_0\nV  Φ₁(A δ)U(k), where Φ₁(A δ) is defined in Φ₁(A, δ; [exp_method]).\n\nHere we allow U to be a sequence of time varying non-deterministic input sets.\n\n\n\n\n\n"
},

{
    "location": "lib/discretize.html#Discretize-1",
    "page": "Discretization",
    "title": "Discretize",
    "category": "section",
    "text": "Pages = [\"discretize.md\"]\nDepth = 3CurrentModule = Reachability.ReachSetsdiscretize\ndiscretize_interpolation\ndiscretize_firstorder\ndiscretize_nobloating"
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
    "location": "publications.html#",
    "page": "Publications",
    "title": "Publications",
    "category": "page",
    "text": ""
},

{
    "location": "publications.html#Publications-1",
    "page": "Publications",
    "title": "Publications",
    "category": "section",
    "text": "Pages = [\"publications.md\"]This page lists publications about the JuliaReach ecosystem and its applications."
},

{
    "location": "publications.html#JuliaReach:-a-Toolbox-for-Set-Based-Reachability-1",
    "page": "Publications",
    "title": "JuliaReach: a Toolbox for Set-Based Reachability",
    "category": "section",
    "text": "JuliaReach: a Toolbox for Set-Based Reachability. Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, Christian Schilling. Accepted in Proceedings of HSCC\'19: 22nd ACM International Conference on Hybrid Systems: Computation and Control (HSCC\'19). Get pdf from arXiv: 1901.10736.In 2019, this conference is part of the Cyber-Physical Systems and Internet-Of-Things Week.Abstract. We present JuliaReach, a toolbox for set-based reachability analysis of dynamical systems. JuliaReach consists of two main packages: Reachability, containing implementations of reachability algorithms for continuous and hybrid systems, and LazySets, a standalone library that implements state-of-the-art algorithms for calculus with convex sets. The library offers both concrete and lazy set representations, where the latter stands for the ability to delay set computations until they are needed. The choice of the programming language Julia and the accompanying documentation of our toolbox allow researchers to easily translate set-based algorithms from mathematics to software in a platform-independent way, while achieving runtime performance that is comparable to statically compiled languages. Combining lazy operations in high dimensions and explicit computations in low dimensions, JuliaReach can be applied to solve complex, large-scale problems.The repeatability evaluation package for this conference tool paper is available at HSCC2019_RE."
},

{
    "location": "publications.html#ARCH-2018-Competition-AFF-Category-Report-1",
    "page": "Publications",
    "title": "ARCH 2018 Competition AFF Category Report",
    "category": "section",
    "text": "ARCH-COMP18 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics. Matthias Althoff, Stanley Bak, Xin Chen, Chuchu Fan,    Marcelo Forets, Goran Frehse, Niklas Kochdumper, Yangge Li, Sayan Mitra, Rajarshi Ray, Christian Schilling and Stefan Schupp (2018) ARCH18. 5th International Workshop on Applied Verification of Continuous and Hybrid Systems, 54: 23–52. doi: 10.29007/73mb. Packages: Reachability.jl.Abstract. This report presents the results of a friendly competition for formal verification of continuous and hybrid systems with linear continuous dynamics. The friendly competition took place as part of the workshop Applied Verification for Continuous and Hybrid Systems (ARCH) in 2018. In its second edition, 9 tools have been applied to solve six different benchmark problems in the category for linear continuous dynamics (in alphabetical order): CORA, CORA/SX, C2E2, FlowStar, HyDRA, Hylaa, Hylaa-Continuous, JuliaReach, SpaceEx, and XSpeed. This report is a snapshot of the current landscape of tools and the types of benchmarks they are particularly suited for. Due to the diversity of problems, we are not ranking tools, yet the presented results probably provide the most complete assessment of tools for the safety verification of continuous and hybrid systems with linear continuous dynamics up to this date.The repeatability evaluation package for JuliaReach is available at ARCH2018_RE.The repeatability evaluation packages of all tools participating in this report is available in the ARCH-COMP gitlab repo."
},

{
    "location": "publications.html#Award-to-JuliaReach-1",
    "page": "Publications",
    "title": "Award to JuliaReach",
    "category": "section",
    "text": "The Best Friendly Competition Result of the 2nd International Competition on Verifying Continuous and Hybrid Systems (ARCH) was given to JuliaReach for the results obtained in ARCH2018_RE for the affine category; see the announcement here:It is our pleasure to announce that Marcelo Forets and Christian Schilling receive the ARCH 2018 Best Friendly Competition Result. They develop the tool JuliaReach, which showed significant improvements for computing reachable sets of linear continuous systems. The award comes with a 500 Euro prize from Bosch. Goran Frehse received the prize from Thomas Heinz of Bosch on their behalf."
},

{
    "location": "publications.html#Reach-Set-Approximation-through-Decomposition-1",
    "page": "Publications",
    "title": "Reach Set Approximation through Decomposition",
    "category": "section",
    "text": "Reach Set Approximation through Decomposition with Low-dimensional Sets and High-dimensional Matrices. Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Frédéric Viry, Andreas Podelski and Christian Schilling (2018) HSCC\'18 Proceedings of the 21st International Conference on Hybrid Systems: Computation and Control: 41–50. See the ACM Digital Library link, or the arXiv: 1801.09526. Packages: LazySets.jl and Reachability.jl. Abstract. Approximating the set of reachable states of a dynamical system is an algorithmic yet mathematically rigorous way to reason about its safety. Although progress has been made in the development of efficient algorithms for affine dynamical systems, available algorithms still lack scalability to ensure their wide adoption in the industrial setting. While modern linear algebra packages are efficient for matrices with tens of thousands of dimensions, set-based image computations are limited to a few hundred. We propose to decompose reach set computations such that set operations are performed in low dimensions, while matrix operations like exponentiation are carried out in the full dimension. Our method is applicable both in dense- and discrete-time settings. For a set of standard benchmarks, it shows a speed-up of up to two orders of magnitude compared to the respective state-of-the art tools, with only modest losses in accuracy. For the dense-time case, we show an experiment with more than 10.000 variables, roughly two orders of magnitude higher than possible with previous approaches.For the models with the SLICOT benchmarks and the repeatability evaluation see ReachabilityBenchmarks."
},

{
    "location": "citations.html#",
    "page": "Citations",
    "title": "Citations",
    "category": "page",
    "text": ""
},

{
    "location": "citations.html#Citations-1",
    "page": "Citations",
    "title": "Citations",
    "category": "section",
    "text": "Pages = [\"citations.md\"]This page lists publications citing packages or papers from the JuliaReach ecosystem."
},

{
    "location": "citations.html#Conference-Proceedings-1",
    "page": "Citations",
    "title": "Conference Proceedings",
    "category": "section",
    "text": "Schupp, Stefan, and Erika Ábrahám. \"Spread the Work: Multi-threaded Safety Analysis for Hybrid Systems.\" International Conference on Software Engineering and Formal Methods. Springer, Cham, 2018.\nBak, Stanley, Hoang-Dung Tran, and Taylor T. Johnson. \"Numerical Verification of Affine Systems with up to a Billion Dimensions.\" arXiv preprint arXiv:1804.01583 (2018). Accepted in Proceedings of HSCC\'19: 22nd ACM International Conference on Hybrid Systems: Computation and Control (HSCC\'19).\nSchupp, Stefan, Justin Winkens, and Erika Ábrahám. \"Context-Dependent Reachability Analysis for Hybrid Systems.\" 2018 IEEE International Conference on Information Reuse and Integration (IRI). IEEE, 2018."
},

{
    "location": "citations.html#Preprints-1",
    "page": "Citations",
    "title": "Preprints",
    "category": "section",
    "text": "Mitchell, Ian M., Jacob Budzis, and Andriy Bolyachevets. \"Invariant, Viability and Discriminating Kernel Under-Approximation via Zonotope Scaling.\" arXiv preprint arXiv:1901.01006 (2019)."
},

{
    "location": "citations.html#Theses-1",
    "page": "Citations",
    "title": "Theses",
    "category": "section",
    "text": "Rocca, Alexandre. Formal methods for modelling and validation of biological models. Diss. Grenoble Alpes, 2018.\nSchilling, Christian. \"Fundamental techniques for the scalable analysis of systems.\" (2018)."
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