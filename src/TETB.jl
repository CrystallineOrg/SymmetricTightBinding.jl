module TETB

# -------------------- Necessary modules for the package --------------------------------- #
using Crystalline
using MPBUtils
# TODO: ↓ probably not necessary to add (`SymmetryBases` and `PhotonicBandConnectivity`) 
#         after your changes in `MPBUtils.jl`
using SymmetryBases
using PhotonicBandConnectivity
#
using Crystalline: AbstractSymmetryVector, irdim, CompositeBandRep_from_indices, translation
using Crystalline: reduce_translation_to_unitrange, constant, free, isapproxin, orbit

# -------------------- Import MPB and Meep ----------------------------------------------- #
using PyCall

const mp = PyNULL()
const mpb = PyNULL()
function __init__()
    # import the mp and mpb libraries
    # https://github.com/JuliaPy/PyCall.jl#using-pycall-from-julia-modules)
    try
        copy!(mp, pyimport("meep"))
        copy!(mpb, pyimport("meep.mpb"))
    catch
        @warn "mpb or meep could not be imported: associated functionality is nonfunctional"
    end
end

export mp, mpb

# ---------------------------------------------------------------------------------------- #

include("SymmetricTightBinding.jl/SymmetricTightBinding.jl")

include("conversion.jl")

include("constrained_nonnegative_expansions.jl")
include("utils.jl")
export find_apolar_modes
export find_auxiliary_modes
export find_bandrep_decompositions
export obtain_symmetry_vectors

# ---------------------------------------------------------------------------------------- #
end # module
