# ===============================
# Test default values and aliases
# ===============================
o = Reachability.Options(:key1 => 1, :key3 => 1, :alias4 => "a")
dict1 = o.dict
dict2 = Dict{Symbol, Any}()

# do not modify dict1 for existing key
Reachability.check_aliases_and_add_default_value!(dict1, dict2, [:key1], 2)
@test haskey(dict1, :key1) && dict1[:key1] == 1
@test haskey(dict2, :key1) && dict2[:key1] == 1

# use default value
Reachability.check_aliases_and_add_default_value!(dict1, dict2, [:key2], 2, true)
@test haskey(dict1, :key2) && dict1[:key2] == 2
@test haskey(dict2, :key2) && dict2[:key2] == 2

# an unused alias has no effect
Reachability.check_aliases_and_add_default_value!(dict1, dict2, [:key3, :alias3], 2)
@test dict1[:key3] == 1
@test !haskey(dict1, :alias3)
@test haskey(dict2, :key3) && dict2[:key3] == 1
@test !haskey(dict2, :alias3)

# check alias detection
Reachability.check_aliases_and_add_default_value!(dict1, dict2, [:key4, :alias4], "b")
@test !haskey(dict1, :key4)
@test haskey(dict1, :alias4) && dict1[:alias4] == "a"
@test haskey(dict2, :key4) && dict2[:key4] == "a"
@test !haskey(dict2, :alias4)

# check alias detection together with default values
Reachability.check_aliases_and_add_default_value!(dict1, dict2, [:key5, :alias5], "b", true)
@test haskey(dict1, :key5) && dict1[:key5] == "b"
@test !haskey(dict1, :alias5)
@test haskey(dict2, :key5) && dict2[:key5] == "b"
@test !haskey(dict2, :alias5)
