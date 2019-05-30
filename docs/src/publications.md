# Publications

```@contents
Pages = ["publications.md"]
```

This page lists publications about the JuliaReach ecosystem and its applications.

## ARCH 2019 Competition NLN Category Report

- **ARCH-COMP19 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics.** Fabian Immler, Matthias Althoff, Luis Benet, Alexandre Chapoutot, Xin Chen, Marcelo Forets, Luca Geretti, Niklas Kochdumper, David P. Sanders and Christian Schilling (2019) ARCH19. 6th International Workshop on Applied Verification of Continuous and Hybrid Systems, vol 61, pages 41--61 [doi: 10.29007/bj1w](https://easychair.org/publications/paper/1gbP). Packages: [Reachability.jl](https://github.com/JuliaReach/Reachability.jl).

**Abstract.** _We present the results of a friendly competition for formal verification of continuous and hybrid systems with nonlinear continuous dynamics. The friendly competition took place as part of the workshop Applied Verification for Continuous and Hybrid Systems (ARCH) in 2019. In this year, 6 tools Ariadne, CORA, DynIbex, Flow*, Isabelle/HOL, and JuliaReach (in alphabetic order) participated. They are applied to solve reachability analysis problems on four benchmark problems, one of them with hybrid dynamics. We do not rank the tools based on the results, but show the current status and discover the potential advantages of different tools._

The repeatability evaluation package for JuliaReach is available at [ARCH2019_RE](https://github.com/JuliaReach/ARCH2019_RE).

The repeatability evaluation packages of all tools participating in this report is available in the [ARCH-COMP gitlab repo](https://gitlab.com/goranf/ARCH-COMP).
 
> [JuliaReach](http://juliareach.org) is a toolbox for reachability computations of dynamical systems, available at http://juliareach.org. For nonlinear reachability we combine functionality from [Taylor-Models.jl](https://github.com/JuliaIntervals/TaylorModels.jl/), [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl) and [TaylorIntegration.jl](https://github.com/PerezHz/TaylorIntegration.jl). First, we compute a non-validated integration using a Taylor model of  order `nT`. The coefficients of that series are polynomials of ordern `Q` in the variables that denote the small variations of the initial conditions. We obtain a time step from the last two coefficients of this time series. In order to validate the  integration step, we compute a second integration using intervals as coefficients of the polynomials in time, and we obtain a bound for the integration using a Lagrange-like remainder. The remainder is used to check the contraction of a Picard iteration. If the combination of the time step and the remainder do not satisfy the contraction, we iteratively enlarge the remainder or possibly shrink the time step. Finally, we evaluate the initial Taylor series with the valid remainder at the time step for which the contraction has been proved, which is also evaluated in the initial box of deviations from the central initial condition, to yield an over-approximation. The approach is (numerically) sound due to rigorous interval bounds in the Taylor approximation.

## ARCH 2019 Competition AFF Category Report

- **ARCH-COMP19 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Marcelo Forets, Goran Frehse, Niklas Kochdumper, Rajarshi Ray, Christian Schilling and Stefan Schupp (2019) ARCH19. 6th International Workshop on Applied Verification of Continuous and Hybrid Systems, vol 61, pages 14--40 [doi: 10.29007/bj1w](https://easychair.org/publications/paper/1gbP). Packages: [Reachability.jl](https://github.com/JuliaReach/Reachability.jl).

**Abstract.** *This report presents the results of a friendly competition for formal verification of continuous and hybrid systems with linear continuous dynamics. The friendly competition took place as part of the workshop Applied Verification for Continuous and Hybrid Systems (ARCH) in 2019. In its third edition, seven tools have been applied to solve six different benchmark problems in the category for linear continuous dynamics (in alphabetical order): CORA, CORA/SX, HyDRA, Hylaa, JuliaReach, SpaceEx, and XSpeed. This report is a snapshot of the current landscape of tools and the types of benchmarks they are particularly suited for. Due to the diversity of problems, we are not ranking tools, yet the presented results provide one of the most complete assessments of tools for the safety verification of continuous and hybrid systems with linear continuous dynamics up to this date.*

The repeatability evaluation package for JuliaReach is available at [ARCH2019_RE](https://github.com/JuliaReach/ARCH2019_RE).

The repeatability evaluation packages of all tools participating in this report is available in the [ARCH-COMP gitlab repo](https://gitlab.com/goranf/ARCH-COMP).

## Reachability analysis of linear hybrid systems via block decomposition

- **Reachability analysis of linear hybrid systems via block decomposition.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, Christian Schilling. [Get pdf from arXiv: 1905.02458](https://arxiv.org/abs/1905.02458).

**Abstract.** *Reachability analysis aims at identifying states reachable by a system within a given time horizon. This task is known to be computationally hard for hybrid systems. One of the main challenges is the handling of discrete transitions, including computation of intersections with invariants and guards. In this paper, we address this problem by proposing a state-space decomposition approach for linear hybrid systems. This approach allows us to perform most operations in low-dimensional state space, which can lead to significant performance improvements.*


## JuliaReach: a Toolbox for Set-Based Reachability

- **JuliaReach: a Toolbox for Set-Based Reachability.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Kostiantyn Potomkin, Christian Schilling. Published in Proceedings of [HSCC'19](http://hscc2019.eecs.umich.edu/): 22nd ACM International Conference on Hybrid Systems: Computation and Control (HSCC'19), see [ACM link here](https://dl.acm.org/citation.cfm?id=3311804). [Get pdf from arXiv: 1901.10736](https://arxiv.org/abs/1901.10736).

In 2019, this conference is part of the [Cyber-Physical Systems and Internet-Of-Things Week](http://www.cpsweek.org/).

**Abstract.** *We present JuliaReach, a toolbox for set-based reachability analysis of
dynamical systems. JuliaReach consists of two main packages: Reachability,
containing implementations of reachability algorithms for continuous and hybrid
systems, and LazySets, a standalone library that implements state-of-the-art
algorithms for calculus with convex sets. The library offers both concrete and
lazy set representations, where the latter stands for the ability to delay set
computations until they are needed. The choice of the programming language
Julia and the accompanying documentation of our toolbox allow researchers to
easily translate set-based algorithms from mathematics to software in a
platform-independent way, while achieving runtime performance that is
comparable to statically compiled languages. Combining lazy operations in high
dimensions and explicit computations in low dimensions, JuliaReach can be
applied to solve complex, large-scale problems.*

The repeatability evaluation package for this conference tool paper is available
at [HSCC2019_RE](https://github.com/JuliaReach/HSCC2019_RE).


## ARCH 2018 Competition AFF Category Report

- **ARCH-COMP18 Category Report: Continuous and Hybrid Systems with Linear Continuous Dynamics.** Matthias Althoff, Stanley Bak, Xin Chen, Chuchu Fan,    Marcelo Forets, Goran Frehse, Niklas Kochdumper, Yangge Li, Sayan Mitra, Rajarshi Ray, Christian Schilling and Stefan Schupp (2018) ARCH18. 5th International Workshop on Applied Verification of Continuous and Hybrid Systems, 54: 23–52. [doi: 10.29007/73mb](http://dx.doi.org/10.29007/73mb). Packages: [Reachability.jl](https://github.com/JuliaReach/Reachability.jl).

**Abstract.** *This report presents the results of a friendly competition for formal verification of continuous and hybrid systems with linear continuous dynamics. The friendly competition took place as part of the workshop Applied Verification for Continuous and Hybrid Systems (ARCH) in 2018. In its second edition, 9 tools have been applied to solve six different benchmark problems in the category for linear continuous dynamics (in alphabetical order): CORA, CORA/SX, C2E2, FlowStar, HyDRA, Hylaa, Hylaa-Continuous, JuliaReach, SpaceEx, and XSpeed. This report is a snapshot of the current landscape of tools and the types of benchmarks they are particularly suited for. Due to the diversity of problems, we are not ranking tools, yet the presented results probably provide the most complete assessment of tools for the safety verification of continuous and hybrid systems with linear continuous dynamics up to this date.*

The repeatability evaluation package for JuliaReach is available at [ARCH2018_RE](https://github.com/JuliaReach/ARCH2018_RE).

The repeatability evaluation packages of all tools participating in this report is available in the [ARCH-COMP gitlab repo](https://gitlab.com/goranf/ARCH-COMP).

### Award to JuliaReach

The *Best Friendly Competition Result* of the [2nd International Competition on Verifying Continuous and Hybrid Systems (ARCH)](https://cps-vo.org/group/ARCH/FriendlyCompetition) was given to JuliaReach for the results obtained in [ARCH2018_RE](https://github.com/JuliaReach/ARCH2018_RE) for the affine category; see the [announcement here](https://cps-vo.org/node/55228):

> It is our pleasure to announce that Marcelo Forets and Christian Schilling receive the ARCH 2018 Best Friendly Competition Result. They develop the tool JuliaReach, which showed significant improvements for computing reachable sets of linear continuous systems. The award comes with a 500 Euro prize from Bosch. Goran Frehse received the prize from Thomas Heinz of Bosch on their behalf.


## Reach Set Approximation through Decomposition

- **Reach Set Approximation through Decomposition with Low-dimensional Sets and High-dimensional Matrices.** Sergiy Bogomolov, Marcelo Forets, Goran Frehse, Frédéric Viry, Andreas Podelski and Christian Schilling (2018) [HSCC'18](https://www.hscc2018.deib.polimi.it/) Proceedings of the 21st International Conference on Hybrid Systems: Computation and Control: 41–50. See the [ACM Digital Library link](http://dx.doi.org/10.1145/3178126.3178128), or the [arXiv: 1801.09526](https://arxiv.org/abs/1801.09526). Packages: [LazySets.jl](https://github.com/JuliaReach/LazySets.jl) and [Reachability.jl](https://github.com/JuliaReach/Reachability.jl). 

**Abstract.** *Approximating the set of reachable states of a dynamical system is an algorithmic yet mathematically rigorous way to reason about its safety. Although progress has been made in the development of efficient algorithms for affine dynamical systems, available algorithms still lack scalability to ensure their wide adoption in the industrial setting. While modern linear algebra packages are efficient for matrices with tens of thousands of dimensions, set-based image computations are limited to a few hundred. We propose to decompose reach set computations such that set operations are performed in low dimensions, while matrix operations like exponentiation are carried out in the full dimension. Our method is applicable both in dense- and discrete-time settings. For a set of standard benchmarks, it shows a speed-up of up to two orders of magnitude compared to the respective state-of-the art tools, with only modest losses in accuracy. For the dense-time case, we show an experiment with more than 10.000 variables, roughly two orders of magnitude higher than possible with previous approaches.*

For the models with the SLICOT benchmarks and the repeatability evaluation see [ReachabilityBenchmarks](https://github.com/JuliaReach/ReachabilityBenchmarks).
