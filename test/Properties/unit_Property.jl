# ====================
# Test Clause creation
# ====================
lc = LinearConstraint([1.; zeros(7)], 0.35)
# vector of atoms
c1 = Clause([lc, LinearConstraint([zeros(4); 1.; zeros(3)], 0.45)])
# single atom
c2 = Clause(lc)

# ======================
# Test Property creation
# ======================
# vector of clauses
p1 = Property([c1, c2])
# single clause
p2 = Property(c1)
# single constraint
p3 = Property(lc)
# single raw constraint
p4 = Property([1.; zeros(7)], 0.35)


# just have one actual test case...
@test true
