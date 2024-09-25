using Pkg
Pkg.activate(@__DIR__)

### necessary packages
using Crystalline
using SymmetryBases, MPBUtils
using PhotonicBandConnectivity
using TETB

### construct the structure under study

R1 = 0.2 #cylinder radius
mat = mp.Medium(epsilon=12)
geometry = [
    mp.Cylinder(radius=R1, center=[0, 0, 0], axis=[0, 0, 1], height=1, material=mat),
    mp.Cylinder(radius=R1, center=[0, 0, 0], axis=[0, 1, 0], height=1, material=mat),
    mp.Cylinder(radius=R1, center=[0, 0, 0], axis=[1, 0, 0], height=1, material=mat),
]

### solve the system
ms = mpb.ModeSolver(
    num_bands=8,
    geometry_lattice=mp.Lattice(basis1=[1, 0, 0], basis2=[0, 1, 0], basis3=[0, 0, 1],
        size=[1, 1, 1]),
    geometry=geometry,
    resolution=16,
)
ms.init_params(p=mp.ALL, reset_fields=true)

### obtain the symmetry vectors of the bands computed above
sg_num = 221
symvecs, topologies = obtain_symmetry_vectors(ms, sg_num)

m = symvecs[1] # pick the 2 lower bands

### obtain additional modes with dimendion `t`
μᴸ = 1
brs = calc_bandreps(sg_num)
idxsᴸs = find_auxiliary_modes(μᴸ, brs)

### compute all possible decomposition into EBRs of m using the additional modes computed
candidatesv = find_apolar_modes(m, idxsᴸs, brs)

#------------------------------------------------------------------------------------------#
# make a TB model out of one of the solutions

## fisrt I need to construct the representations of the operations of the SG. At least the 
## generators
nᵀ⁺ᴸ = candidatesv[1][1][1]

# TODO: I did this manually, can I extract it from nᵀ⁺ᴸ?
gen = generators(sg_num)
wp = wyckoffs(221)[end-2]
siteg = sitegroup(sg_num, wp)
cosetsg = cosets(siteg)

siterep = nᵀ⁺ᴸ.siteir

n_wp = wp.mult
dim_rep = size(siterep.matrices[1][1]) # TODO: maybe some dim of the siteir?

wps = orbit(group(nᵀ⁺ᴸ))
siteir = nᵀ⁺ᴸ.siteir
siteir_dim = Crystalline.irdim(siteir)
##

using Crystalline: irdim, constant

# we do not include the (usually redundant) exponential phases below
# TODO: check & write doc string describing what this does
function sgrep_induced_by_siteir_generators(br::NewBandRep{D}) where D
    siteir = br.siteir
    siteir_dim = irdim(siteir)
    siteg = group(siteir)
    wps = orbit(siteg)
    mult = length(wps)
    gens = generators(num(siteg))
    
    ρs = [BlockArray{ComplexF64}(
            zeros(ComplexF64, siteir_dim*mult, siteir_dim*mult),
            fill(siteir_dim, mult), fill(siteir_dim, mult)) for _ in eachindex(gens)]
    for (n, g) in enumerate(gens)
        ρ = ρs[n]
        for (α, (gₐ, qₐ)) in enumerate(zip(cosets(siteg), wps))
            check = false
            for (β, (gᵦ, qᵦ)) in enumerate(zip(cosets(siteg), wps))
                tᵦₐ = constant(g*parent(qₐ) - parent(qᵦ)) # ignore free parts
                # compute h = gᵦ⁻¹ tᵦₐ⁻¹ g gₐ
                h = compose(compose(compose(inv(gᵦ), SymOperation(-tᵦₐ), false), g, false), gₐ, false)
                idx_h = findfirst(==(h), siteg)
                if !isnothing(idx_h) # h ∈ siteg and qₐ and qᵦ are connected by g
                    ρ[Block(α, β)] .= siteir.matrices[idx_h]
                    check = true
                    break
                end
            end
            check || error("failed to find any nonzero block")
        end
    end

    return gens .=> ρs
end
br = brs[1]
sgrep_induced_by_siteir_generators(br)[1][2]

# for i in gen
#     D = zeros((n_wp, dim_rep))


# end

#------------------------------------------------------------------------------------------#

