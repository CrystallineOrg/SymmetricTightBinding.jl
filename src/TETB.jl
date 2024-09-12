module TETB

# ---------------------------------------------------------------------------------------- #

using Crystalline
using MPBUtils
using SymmetryBases
using PhotonicBandConnectivity
const PBC = PhotonicBandConnectivity

# ---------------------------------------------------------------------------------------- #

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

include("read_utils.jl")
include("constrained_nonnegative_expansions.jl")
include("utils.jl")

export find_all_band_representations, find_auxiliary_modes, obtain_symmetry_vectors, find_physical_band_representations

# ---------------------------------------------------------------------------------------- #
end # module
