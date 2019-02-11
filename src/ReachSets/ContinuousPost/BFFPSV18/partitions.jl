"""
    inout_map_property(prop::IntersectionProperty,
                       partition::AbstractVector{<:AbstractVector{Int}},
                       blocks::AbstractVector{Int},
                       n::Int
                      )::IntersectionProperty

Map an `IntersectionProperty` to the dimensions of analyzed blocks.

### Input

- `prop`      -- property
- `partition` -- block partition; elements are start and end indices of a block
- `blocks`    -- list of all output block indices in the partition
- `n`         -- total number of input dimensions

### Output

A new property of reduced dimension.

### Notes

If the dimension is not reduced, we keep the original set.
Otherwise, the dimension reduction is achieved with a `LinearMap`.
"""
function inout_map_property(prop::IntersectionProperty{N},
                            partition::AbstractVector{<:AbstractVector{Int}},
                            blocks::AbstractVector{Int},
                            n::Int
                           )::IntersectionProperty{N} where {N<:Real}
    @assert dim(prop.bad) == n "the property has dimension $(dim(prop.bad)) but should have dimension $n"
    proj = projection_map(partition, blocks)
    if length(proj) == n
        # no change in the dimension, copy the old property (keep the set)
        return IntersectionProperty(prop.bad)
    else
        M = sparse(proj, proj, ones(N, length(proj)), n, n)
        return IntersectionProperty(M * prop.bad)
    end
end

"""
    inout_map_property(prop::LinearConstraintProperty{N},
                       partition::AbstractVector{<:AbstractVector{Int}},
                       blocks::AbstractVector{Int},
                       n::Int
                      )::LinearConstraintProperty{N} where {N<:Real}

Map a `LinearConstraintProperty` to the dimensions of analyzed blocks.

### Input

- `prop`      -- property
- `partition` -- block partition; elements are start and end indices of a block
- `blocks`    -- list of all output block indices in the partition
- `n`         -- total number of input dimensions

### Output

A new property of reduced dimension.
"""
function inout_map_property(prop::LinearConstraintProperty{N},
                            partition::AbstractVector{<:AbstractVector{Int}},
                            blocks::AbstractVector{Int},
                            n::Int
                           )::LinearConstraintProperty{N} where {N<:Real}
    # sanity check: do not project away non-zero dimensions
    function check_projection(a, proj)
        p = 1
        for i in 1:length(a)
            if p <= length(proj) && i == proj[p]
                # dimension is not projected away
                p += 1
            elseif a[i] != 0
                # dimension is projected away; entry is non-zero
                return false
            end
        end
        return true
    end

    @assert dim(prop.clauses[1].atoms[1]) == n "the property has dimension $(dim(prop.clauses[1].atoms[1])) but should have dimension $n"

    proj = projection_map(partition, blocks)

    # create modified property
    clauses = Vector{Clause{N}}(undef, length(prop.clauses))
    for (ic, c) in enumerate(prop.clauses)
        atoms = Vector{LinearConstraint{N}}(undef, length(c.atoms))
        for (ia, atom) in enumerate(c.atoms)
            @assert check_projection(atom.a, proj) "blocks incompatible with property"
            atoms[ia] = LinearConstraint{N}(atom.a[proj], atom.b)
        end
        clauses[ic] = Clause(atoms)
    end
    return LinearConstraintProperty(clauses)
end

"""
    projection_map(partition::AbstractVector{<:AbstractVector{Int}},
                   blocks::AbstractVector{Int}
                  )::Vector{Int}

Create a sorted list of all dimensions in `blocks`.

### Input

- `partition` -- block partition; elements are start and end indices of a block
- `blocks`    -- list of all output block indices in the partition

### Output

A sorted list containing all output dimensions.
"""
function projection_map(partition::AbstractVector{<:AbstractVector{Int}},
                        blocks::AbstractVector{Int}
                       )::Vector{Int}
    proj = Vector{Int}()
    for bi in blocks
        for i in partition[bi]
            push!(proj, i)
        end
    end
    return proj
end

"""
    inout_map_property(prop::SubsetProperty,
                       partition::AbstractVector{<:AbstractVector{Int}},
                       blocks::AbstractVector{Int},
                       n::Int
                      )::SubsetProperty

Map an `SubsetProperty` to the dimensions of analyzed blocks.

### Input

- `prop`      -- property
- `partition` -- block partition; elements are start and end indices of a block
- `blocks`    -- list of all output block indices in the partition
- `n`         -- total number of input dimensions

### Output

A new property of reduced dimension.

### Notes

If the dimension is not reduced, we keep the original set.
Otherwise, the dimension reduction is achieved with a `LinearMap`.
"""
function inout_map_property(prop::SubsetProperty{N},
                            partition::AbstractVector{<:AbstractVector{Int}},
                            blocks::AbstractVector{Int},
                            n::Int
                           )::SubsetProperty{N} where {N<:Real}
    @assert dim(prop.safe) == n "the property has dimension $(dim(prop.safe)) but should have dimension $n"
    proj = projection_map(partition, blocks)
    if length(proj) == n
        # no change in the dimension, copy the old property (keep the set)
        return IntersectionProperty(prop.safe)
    else
        M = sparse(proj, proj, ones(N, length(proj)), n, n)
        return IntersectionProperty(M * prop.safe)
    end
end
