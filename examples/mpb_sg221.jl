using Pkg
Pkg.activate(@__DIR__)

### necessary packages
using Crystalline
using SymmetryBases, MPBUtils
using PhotonicBandConnectivity
using TETB
using BlockArrays

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

sgrep = sgrep_induced_by_siteir_generators(nᵀ⁺ᴸ) # representation of the SG generators in the
# basis defined by nᵀ⁺ᴸ

#------------------------------------------------------------------------------------------#

