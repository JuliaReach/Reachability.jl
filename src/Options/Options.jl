__precompile__()
"""
Module to handle options of reachability algorithms.
"""
module Options

    #=============================================================
    Struct with the dictionary of options and basic functionality
    =============================================================#
    include("dictionary.jl")

    #=============================================
    Validation methods for dictionaries of options
    ==============================================#
    include("validation.jl")
end
