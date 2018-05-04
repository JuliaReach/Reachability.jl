## Introduction

The SLICOT (for Subroutine Library in Systems and Control Theory) benchmarks are a set of mathematical models which
reflect “real world” applications [1], [2], [3]. The models are publicly available,
see the [Benchmark Examples for Model Reduction](http://slicot.org/20-site/126-benchmark-examples-for-model-reduction)
webpage for further details and to download the models.

## Results

In the project [ReachabilityBenchmarks.jl](https://github.com/JuliaReach/ReachabilityBenchmarks/) we have
written Julia scripts to run the models `motor` (8 variables), `building` (48 variables),
`cdplayer`,
Some of the original models are differential algebraic equations (DAEs), in which case we only kept the ODE part,
i.e., the coefficient matrices A and B, which is consistent with related literature on reach set approximation.

## References

[1] P. Benner, V. Mehrmann, V. Sima, S. Van Huffel, and A. Varga, “SLICOT – a subroutine library in systems and control theory,” in Applied and computational control, signals, and circuits. Springer, 1999, pp. 499–539.

[2] Y. Chahlaoui and P. Van Dooren, “Benchmark examples for model reduction of linear time-invariant dynamical systems,” in Dimension Reduction of Large-Scale Systems. Springer, 2005, pp. 379–392.

[3] H. Tran, L. V. Nguyen, and T. T. Johnson, “Large-scale linear systems from order-reduction,” in ARCH, vol. 43. EasyChair, 2016, pp. 60–67.
