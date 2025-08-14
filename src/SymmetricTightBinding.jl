module SymmetricTightBinding

# --- Dependencies ----------------------------------------------------------------------- #

using LinearAlgebra
using Crystalline
using Crystalline: AbstractSymmetryVector, irdim, CompositeBandRep_from_indices, translation
using Crystalline: reduce_translation_to_unitrange, constant, free, isapproxin, orbit
using BlockArrays
using RowEchelon: rref, rref!           # for `poormans_sparsification`

# --- Constants -------------------------------------------------------------------------- #

const VEC_CMP_ATOL = 1e-11 # for `isapprox` comparison of RVecs / KVecs
const NULLSPACE_ATOL_DEFAULT = 1e-5
const SPARSIFICATION_ATOL_DEFAULT = 1e-10
const PRUNE_ATOL_DEFAULT = SPARSIFICATION_ATOL_DEFAULT
const ZASSENHAUS_ATOL_DEFAULT = NULLSPACE_ATOL_DEFAULT

# --- Code loading ----------------------------------------------------------------------- #

include("types.jl")
export HoppingOrbit
export TightBindingElementString
export TightBindingBlock
export TightBindingModel
export ParameterizedTightBindingModel
include("show.jl")
include("site_representations.jl")
export sgrep_induced_by_siteir
include("tightbinding.jl")
export obtain_symmetry_related_hoppings
export tb_hamiltonian
include("zassenhaus.jl")
include("timereversal.jl")
include("hermiticity.jl")
include("utils.jl")
export pin_free!
include("symmetry_analysis.jl")
export symmetry_eigenvalues
include("spectrum.jl")
export spectrum
include("gradients.jl")
export gradient_wrt_hopping
export energy_gradient_wrt_hopping
include("symmetry_breaking.jl")
export subduced_complement

# --- Re-exports ------------------------------------------------------------------------- #

export collect_compatible, collect_irrep_annotations # extended functions from Crystalline

# --- Function defs. & exports for extensions -------------------------------------------- #

function fit end # for Optim.jl extension
export fit

# ---------------------------------------------------------------------------------------- #
end # module