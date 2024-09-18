using Pkg
Pkg.activate(@__DIR__)

using Crystalline, MPBUtils
using PhotonicBandConnectivity, SymmetryBases
using TETB

### construct the structure under srudy

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
    resolution=16,
)
ms.init_params(p=mp.ALL, reset_fields=true)

### obtain the symmetry vectors of the bands under study
sg_num = 221
symvecs, topologies = obtain_symmetry_vectors(ms, sg_num)

vᵀ = symvecs[1] # pick the 2 lower bands

t = 1
brs = calc_bandreps(sg_num)
d = stack(brs)[end, :]

long_modes = find_auxiliary_modes(t, d, brs)

all_band_repre = find_all_band_representations(vᵀ, long_modes, d, brs)