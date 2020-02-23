"""
    inout_map_property(𝑃::Predicate,
                       partition::AbstractVector{<:AbstractVector{Int}},
                       blocks::AbstractVector{Int},
                       n::Int
                      )

Map a property to the dimensions of analyzed blocks.

### Input

- `𝑃`         -- property
- `partition` -- block partition; elements are start and end indices of a block
- `blocks`    -- list of all output block indices in the partition
- `n`         -- total number of input dimensions

### Output

A new property of reduced dimension, if necessary.

### Notes

The return type depends on the ambient dimension (`n`), the block indices in the
partition and on the type of property:

- The property `𝑃` is returned unchanged whenever `n` matches the dimensions
  in `blocks`.
- If the property is a `HalfSpace` (resp. `HPolyhedron`), the function `project`
  from `LazySets` returns a new property with a `HalfSpace` (resp. `HPolyhedron`)
  in the reduced dimensions, according to `blocks`.
- Otherwise, the dimensional reduction is implemented via a (lazy) `LinearMap`.
"""
function inout_map_property(𝑃::Predicate,
                            partition::AbstractVector{<:AbstractVector{Int}},
                            blocks::AbstractVector{Int},
                            n::Int
                           )

    # create a sorted list of all dimensions in `blocks` => available variables for projection
    proj = vcat(partition[blocks]...)

    # no change in the dimension => return the original property
    length(proj) == n && return 𝑃

    return project(𝑃, proj)
end
