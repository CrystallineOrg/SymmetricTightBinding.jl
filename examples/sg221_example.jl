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
    resolution=32,
)
ms.init_params(p=mp.ALL, reset_fields=true)

### obtain the symmetry vectors of the bands under study
sg_num = 221
band_summaries = obtain_symmetry_vectors(ms, sg_num)

vᵀ = band_summaries[1] # pick the 2 lower bands

t = 2
brs = bandreps(sg_num)
d = matrix(brs)[end, :]

long_modes = find_auxiliary_modes(t, d, brs)

all_band_repre = find_all_band_representations(vᵀ, long_modes, d, brs, sg_num)

# TODO: it needs to be generalized to multiple band representations
for i in 1:length(all_band_repre.long_modes)
    nᴸ = [brs[j...] for j in all_band_repre.long_modes[i]]
    print("Solutions using the auxiliary mode: ")

    for j in nᴸ[1:end-1]
        print(j.label, " at ", j.wyckpos, " ⊕ ")
    end

    println(nᴸ[end].label, " at ", nᴸ[end].wyckpos)

    nᵀ⁺ᴸ = [[brs[m...] for k in j for m in k] for j in all_band_repre.solutions[i]]
    count = 1
    for j in nᵀ⁺ᴸ
        print("   ↪Solution #$count: ")
        for k in j[1:end-1]
            print(k.label, " at ", k.wyckpos, " ⊕ ")
        end
        println(j[end].label, " at ", j[end].wyckpos)
        count += 1
    end
end

println("nᵀ⁺ᴸ", " = ", nᵀ⁺ᴸ.label, " at ", nᵀ⁺ᴸ.wyckpos, "; nᴸ", " = ", nᴸ.label, " at ",
    nᴸ.wyckpos, "; Are they physical? ", phys)

phys_band_repre = find_physical_band_representations(vᵀ, long_modes, d, brs, sg_num)

nᵀ⁺ᴸ = brs[phys_band_repre[1][1]...]
nᴸ = brs[phys_band_repre[1][2]...]

println("nᵀ⁺ᴸ", " = ", nᵀ⁺ᴸ.label, " at ", nᵀ⁺ᴸ.wyckpos, "; nᴸ", " = ", nᴸ.label, " at ",
    nᴸ.wyckpos)