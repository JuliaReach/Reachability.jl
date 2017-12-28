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

# check that x1 + x2 <= 2 doesn't hold
s = solve(ContinuousSystem(A, X0), :T=>0.1, :mode=>"check", :blocks=>[1],
          :property=>LinearConstraintProperty([1., 1.], 2.))
@test s.violation == 1

# check that x1 - x2 <= 2 holds
s = solve(ContinuousSystem(A, X0), :T=>0.1, :mode=>"check", :blocks=>[1],
          :property=>LinearConstraintProperty([1., -1.], 2.))
@test s.violation == -1
