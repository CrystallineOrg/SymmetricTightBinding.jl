module TETB

# -------------------- Necessary modules for the package ---------------------------------- #
using LinearAlgebra: nullspace, norm
using Crystalline
using MPBUtils
using SymmetryBases
using PhotonicBandConnectivity
using Crystalline: AbstractSymmetryVector, irdim, CompositeBandRep_from_indices, translation
using Crystalline: reduce_translation_to_unitrange, constant, free, isapproxin, orbit
using BlockArrays
using RowEchelon: rref, rref!           # for `poormans_sparsification`
using GLMakie
# -------------------- Predefined constant used ------------------------------ #
const NULLSPACE_ATOL_DEFAULT = 1e-5
const SPARSIFICATION_ATOL_DEFAULT = 1e-10
const PRUNE_ATOL_DEFAULT = SPARSIFICATION_ATOL_DEFAULT
const ZASSENHAUS_ATOL_DEFAULT = NULLSPACE_ATOL_DEFAULT
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

# ---------------------------------------------------------------------------- #

include("types.jl")
export TightBindingCandidateSet
export HoppingOrbit
export TightBindingElementString, TightBindingBlock
include("conversion.jl")
include("show.jl")
include("constrained_nonnegative_expansions.jl")
include("utils.jl")
export find_apolar_modes
export find_auxiliary_modes
include("site_representations.jl")
export find_bandrep_decompositions
export obtain_symmetry_vectors
export sgrep_induced_by_siteir_generators
include("tightbinding.jl")
export obtain_symmetry_related_hoppings
export tb_hamiltonian
include("plotting_utils.jl")
export hop_plot
include("zassenhaus.jl")
include("time_reversal.jl")

# ---------------------------------------------------------------------------- #
end # module
