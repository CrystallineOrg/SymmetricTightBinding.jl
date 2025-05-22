module TETB

# --- Dependencies ----------------------------------------------------------------------- #

using Crystalline
using MPBUtils
# TODO: ↓ probably unnecessary dependency
using PhotonicBandConnectivity
using Crystalline: AbstractSymmetryVector, irdim, CompositeBandRep_from_indices, translation
using Crystalline: reduce_translation_to_unitrange, constant, free, isapproxin, orbit
using PythonCall: pynew, pycopy!, pyimport, Py, pyconvert
using Reexport
@reexport using SymmetricTightBinding

# --- PythonCall init -------------------------------------------------------------------- #

using PythonCall

const mp = pynew()
const mpb = pynew()
function __init__()
    try # import the mp and mpb libraries
        pycopy!(mp, pyimport("meep"))
        pycopy!(mpb, pyimport("meep.mpb"))
    catch
        @warn "mpb or meep could not be imported: associated functionality is nonfunctional"
    end
end

export mp, mpb

# --- Code loading ----------------------------------------------------------------------- #

include("conversion.jl")
include("constrained_nonnegative_expansions.jl")
include("utils.jl")
export find_apolar_modes
export find_auxiliary_modes
export find_bandrep_decompositions
export obtain_symmetry_vectors

# ---------------------------------------------------------------------------- #
end # module
