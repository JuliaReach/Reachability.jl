"""
Abstract supertype of properties that can be checked.

Every concrete subtype should provide the following functions:
  - `inout_map_property(::Property, ::AbstractVector{<:AbstractVector{Int}},
                        ::AbstractVector{Int}, ::Int)::Property`
  - `check_property(::LazySet, ::Property)::Bool`
"""
abstract type Property end

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
