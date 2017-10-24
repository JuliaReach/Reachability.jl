export Clause, Property, check_property

"""
Type that represents a disjunction of linear constraints.

FIELDS :

- ``atoms`` -- a vector of linear constraints
"""
struct Clause
    atoms::Vector{LinearConstraint}

    Clause(atoms::Vector{LinearConstraint}) = new(atoms)

    Clause(atom::LinearConstraint) = new([atom])
end

"""
Type that represents a property that can be checked.
A property represents a Boolean function of linear constraints.
The function is given in CNF (conjunctive normal form), i.e., a conjunction of disjunctions of linear constraints.

FIELDS :

- ``clauses`` -- a vector of `Clause` objects
"""
struct Property
    clauses::Vector{Clause}

    Property(clauses::Vector{Clause}) = new(clauses)

    Property(clause::Clause) = new([clause])

    Property(linConst::LinearConstraint) = new([Clause(linConst)])

    Property(linComb::Vector{Float64}, bound::Float64) = new([Clause(LinearConstraint(linComb, bound))])
end

"""
    check_property(sf, prop)

Checks whether a support function satisfies a property.

INPUTS :
- ``sf`` -- the support function
- ``prop`` -- the property to check
"""
@inline function check_property(sf::LazySet, prop::Property)::Bool
    @assert (length(prop.clauses) > 0) "Empty properties are not allowed."
    for clause in prop.clauses
        @assert (length(clause.atoms) > 0) "Empty clauses are not allowed."
        clause_sat = false
        for atom in clause.atoms
            if ρ(atom.a, sf) < atom.b
                clause_sat = true
                break
            end
        end
        if !clause_sat
            return false
        end
    end
    return true
end
