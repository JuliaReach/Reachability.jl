# ===============================
# Test minimal (default) options
# ===============================
using Reachability
A = randn(4, 4); X0 = BallInf(ones(4), 0.1)
s = solve(ContinuousSystem(A, X0), :T=>0.1);

# ===============================
# Test projection
# ===============================
s = solve(ContinuousSystem(A, X0), :T=>0.1, :blocks=>[1,2], :apply_projection=>false);
ps = project(s);

# ===============
# Test check mode
# ===============

A = [ 0.0509836  0.168159  0.95246   0.33644
      0.42377    0.67972   0.129232  0.126662
      0.518654   0.981313  0.489854  0.588326
      0.38318    0.616014  0.518412  0.778765]

# check that x1 + x2 <= 2 doesn't hold
s = solve(ContinuousSystem(A, X0), :T=>0.1, :mode=>"check", :blocks=>[1],
          :property=>LinearConstraintProperty([1., 1.], 2.))
@test s.violation == 1

# check that x1 - x2 <= 2 holds
s = solve(ContinuousSystem(A, X0), :T=>0.1, :mode=>"check", :blocks=>[1],
          :property=>LinearConstraintProperty([1., -1.], 2.))
@test s.violation == -1

# ===============================
# Test reachability options
# ===============================
s = solve(ContinuousSystem(A, X0), :T=>0.1, :lazy_sih=>true);
s = solve(ContinuousSystem(A, X0), :T=>0.1, :lazy_sih=>false);
