# ======================
# Test property creation
# ======================
X1 = LinearConstraint([1.; zeros(7)], 0.35)
X2 = LinearConstraint([zeros(4); 1.; zeros(3)], 0.45)

p1 = SubsetProperty(X1)
p2 = IntersectionProperty(X2)
p3 = Conjunction([p1, p2])
p4 = Disjunction([p1, p2])
p5 = Conjunction([p3, Disjunction([p1, p4]), p2])

# just have one actual test case...
@test true
