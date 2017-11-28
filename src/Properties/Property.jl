export Clause, Property, check_property

"""
Type that represents a disjunction of linear constraints.

FIELDS :

- ``atoms`` -- a vector of linear constraints
"""
struct Clause{N<:Real}
    atoms::Vector{LinearConstraint{N}}
end
# constructor from a single atom
Clause(atom::LinearConstraint{N}) where {N<:Real} = Clause{N}([atom])

"""
Type that represents a property that can be checked.
A property represents a Boolean function of linear constraints.
The function is given in CNF (conjunctive normal form), i.e., a conjunction of disjunctions of linear constraints.

FIELDS :

- ``clauses`` -- a vector of `Clause` objects
"""
struct Property{N<:Real}
    clauses::Vector{Clause{N}}
end
# constructor from a single clause
Property(clause::Clause{N}) where {N<:Real} = Property{N}([clause])
# constructor from a single constraint
Property(linConst::LinearConstraint{N}) where {N<:Real} =
    Property{N}([Clause{N}(linConst)])
# constructor from a single constraint in raw form ax+b
Property(linComb::Vector{N}, bound::N) where {N<:Real} =
    Property{N}([Clause{N}(LinearConstraint{N}(linComb, bound))])

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
