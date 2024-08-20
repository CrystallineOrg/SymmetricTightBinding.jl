using Pkg
Pkg.activate(".")

using Crystalline, MPBUtils
using PyCall, Pkg
using PhotonicBandConnectivity, SymmetryBases

include("extra_functions.jl")

const PBC = PhotonicBandConnectivity # just to shorten things up

mp = pyimport("meep")
mpb = pyimport("meep.mpb")

### Compute the symmetry vector for an especific structure

R1 = 0.2 #cylinder radius
mat = mp.Medium(epsilon=12)
geometry = [
    mp.Cylinder(radius=R1, center=[0, 0, 0], axis=[0, 0, 1], height=1, material=mat),
    mp.Cylinder(radius=R1, center=[0, 0, 0], axis=[0, 1, 0], height=1, material=mat),
    mp.Cylinder(radius=R1, center=[0, 0, 0], axis=[1, 0, 0], height=1, material=mat),
]
ms = mpb.ModeSolver(
    num_bands=8,
    geometry_lattice=mp.Lattice(basis1=[1, 0, 0], basis2=[0, 1, 0], basis3=[0, 0, 1],
        size=[1, 1, 1]),
    geometry=geometry,
    resolution=32,
)
ms.init_params(p=mp.ALL, reset_fields=true)

sg_num = 221                                        # P-4 (Z₂×Z₂ symmetry indicator group)
brs = bandreps(sg_num)                             # elementary band representations
lgs = littlegroups(sg_num)                         # little groups
filter!(((klab, _),) -> klab ∈ klabels(brs), lgs) # restrict to k-points in `brs`
map!(lg -> primitivize(lg, false), values(lgs))   # convert to primitive setting (without reducing translations)
lgirsd = pick_lgirreps(lgs; timereversal=true);    # small irreps associated with `lgs`

symeigsd = Dict{String,Vector{Vector{ComplexF64}}}()
for (klab, lg) in lgs
    kv = mp.Vector3(position(lg)()...)
    ms.solve_kpoint(kv)

    symeigsd[klab] = [Vector{ComplexF64}(undef, length(lg)) for n in 1:ms.num_bands]
    for (i, gᵢ) in enumerate(lg)
        W = mp.Matrix(eachcol(rotation(gᵢ))...) # decompose gᵢ = {W|w}
        w = mp.Vector3(translation(gᵢ)...)
        symeigs = ms.compute_symmetries(W, w)   # compute ⟨Eₙₖ|gᵢDₙₖ⟩ for all bands
        setindex!.(symeigsd[klab], symeigs, i)  # update container of symmetry eigenvalues
    end
end

# --- fix singular photonic symmetry content at Γ, ω=0 --- # TODO: Maybe we can skip this step or try to change the condition on Γ
fixup_gamma_symmetry!(symeigsd, lgs)

# --- analyze connectivity and topology of symmetry data ---
summaries = analyze_symmetry_data(symeigsd, lgirsd, brs)
vᵀ = summaries[1]
vᵀ´ = prune_klab_irreps_vecs(vᵀ, "Γ")

### Now we want to find a possible EBR decomposition of that vector
B = matrix(brs) # matrix of EBRs with their dimensions in the last row
A = B[1:end-1, :]
d = B[end, :]

#### We need also the same matrice but without Γ
brs´ = prune_klab_irreps_brs(brs, "Γ")
B´ = matrix(brs´) # matrix of EBRs with their dimensions in the last row
A´ = B´[1:end-1, :]

#### find all posible longitudinal modes od dimension `t`
t = 1
long_cand = PBC.filling_symmetry_constrained_expansions(t, Int[], d, brs, Int[]); # search for all long. modes of dimension
# `t``
vᴸ = sum(brs[long_cand[1]])
vᴸ´ = sum(brs´[long_cand[1]])

vᵀ⁺ᴸ = vᵀ.n + vᴸ
vᵀ⁺ᴸ´ = vᵀ´.n + vᴸ´

#### solve the actual problem

μᵀ⁺ᴸ = vᵀ⁺ᴸ[end]

brs_cand = PBC.filling_symmetry_constrained_expansions(μᵀ⁺ᴸ, vᵀ⁺ᴸ´, d, brs´, collect(1:size(B´, 1)))

aᵀ⁺ᴸ = sum(brs[brs_cand[1]])
vᴸ = sum(brs[long_cand[1]])
