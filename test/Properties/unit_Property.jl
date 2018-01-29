# ====================
# Test Clause creation
# ====================
lc = LinearConstraint([1.; zeros(7)], 0.35)
# vector of atoms
c1 = Clause([lc, LinearConstraint([zeros(4); 1.; zeros(3)], 0.45)])
# single atom
c2 = Clause(lc)

# ======================
# Test LinearConstraintProperty creation
# ======================
# vector of clauses
p1 = LinearConstraintProperty([c1, c2])
# single clause
p2 = LinearConstraintProperty(c1)
# single constraint
p3 = LinearConstraintProperty(lc)
# single raw constraint
p4 = LinearConstraintProperty([1.; zeros(7)], 0.35)


# just have one actual test case...
@test true
