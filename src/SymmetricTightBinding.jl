module SymmetricTightBinding

# --- Dependencies ----------------------------------------------------------------------- #

using LinearAlgebra
using Crystalline
using Crystalline: AbstractSymmetryVector, irdim, CompositeBandRep_from_indices, translation
using Crystalline: reduce_translation_to_unitrange, constant, free, isapproxin, orbit
using BlockArrays
using RowEchelon: rref, rref!           # for `poormans_sparsification`
using StaticArrays: SVector

# --- Constants -------------------------------------------------------------------------- #

const VEC_CMP_ATOL = 1e-11 # for `isapprox` comparison of RVecs / KVecs
const NULLSPACE_ATOL_DEFAULT = 1e-5
const SPARSIFICATION_ATOL_DEFAULT = 1e-10
const PRUNE_ATOL_DEFAULT = SPARSIFICATION_ATOL_DEFAULT
const ZASSENHAUS_ATOL_DEFAULT = NULLSPACE_ATOL_DEFAULT

# --- Code loading ----------------------------------------------------------------------- #

include("types.jl")
export HoppingOrbit
export TightBindingBlock
export Hermiticity, HERMITIAN, ANTIHERMITIAN, NONHERMITIAN
export hermiticity
export TightBindingModel
export ParameterizedTightBindingModel
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
export spectrum, spectrum_single_k
include("gradients.jl")
export gradient_wrt_hopping
export TightBindingModelHoppingGradient
export energy_gradient_wrt_hopping
export gradient_wrt_momentum
export TightBindingModelMomentumGradient
export energy_gradient_wrt_momentum
include("caching.jl")
export TightBindingCache
include("berry.jl")
export berrycurvature
export chern
export chern_fukui
include("dos.jl")
export densityofstates
include("symmetry_breaking.jl")
export subduced_complement
include("design.jl")
export uniform_kmesh
export IrrepTarget, locate_multiplet
export IrrepIsolationObjective
export isolation_report
include("show.jl")

# --- Re-exports ------------------------------------------------------------------------- #

export collect_compatible, collect_irrep_annotations # extended functions from Crystalline

# --- Function defs. & exports for extensions -------------------------------------------- #

# fitting functionality; all implemented in the Optim.jl extension
function fit end
function multistart_fit end     # ┐ stubs, overloaded in the extension, so that the fitting
function make_fit_objective end # │ machinery can be reused from dependent packages (e.g.,
function spectralmoments end    # ┘ PhotonicTightBinding.jl) without `Base.get_extension`
export fit, multistart_fit, make_fit_objective

# objective-driven design (cf. `src/design.jl`); driver implemented in the Optim.jl extension
function isolate_irrep end
export isolate_irrep

# ---------------------------------------------------------------------------------------- #
end # module