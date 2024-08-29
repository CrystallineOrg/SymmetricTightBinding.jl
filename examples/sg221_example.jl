using Pkg; Pkg.activate(@__DIR__)

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
    resolution=32,
)
ms.init_params(p=mp.ALL, reset_fields=true)

### obtain the symmetry vectors of the bands under study
sg_num = 221
band_summaries = obtain_symmetry_vectors(ms, sg_num)

vᵀ = band_summaries[1] # pick the 2 lower bands

t = 1
brs = bandreps(sg_num)
d = matrix(brs)[end, :]

long_modes = find_auxiliary_modes(t, d, brs)

# vᵀ⁺ᴸ´ = vᵀ´.n + nᴸ´
# μᵀ⁺ᴸ = vᵀ⁺ᴸ´[end]

band_repre = find_all_band_representations(vᵀ, long_modes, d, brs)

nᵀ⁺ᴸ = brs[band_repre[1][1][1]][1]

    nᵀ⁺ᴸ - vᵀ.n

println((nᵀ⁺ᴸ.label, nᴸ.label))