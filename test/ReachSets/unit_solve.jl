using Reachability, Plots
A = randn(4, 4); X0 = BallInf(ones(4), 0.1)
s = solve(ContinuousSystem(A, X0), :T=>0.1);
