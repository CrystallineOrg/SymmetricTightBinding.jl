using Pkg
Pkg.activate(@__DIR__)

### necessary packages
using TETB
using Crystalline
using TETB.PythonCall: pylist, pyconvert
using Brillouin

### construct the structure under study
sgnum = 221 # space group number
Rs = directbasis(sgnum, Val(3))

R1 = 0.2 #cylinder radius
mat = mp.Medium(; epsilon = 12)
geometry = map([[0, 0, 1], [0, 1, 0], [1, 0, 0]]) do axis
    mp.Cylinder(; radius = R1, center = [0, 0, 0], axis = axis, height = 1, material = mat)
end

### solve the system
ms = mpb.ModeSolver(;
    num_bands = 8,
    geometry_lattice = mp.Lattice(; basis1 = Rs[1], basis2 = Rs[2], basis3 = Rs[3]),
    geometry = pylist(geometry),
    resolution = 16,
)
ms.init_params(; p = mp.ALL, reset_fields = true)

### obtain the symmetry vectors of the bands computed above
symvecs = obtain_symmetry_vectors(ms, sgnum)

m = symvecs[1] # pick the 2 lower bands

### obtain an EBR decomposition with at least one additional band
μᴸ_min = 1
brs = calc_bandreps(sgnum)
candidatesv = find_bandrep_decompositions(m, brs; μᴸ_min)

##-----------------------------------------------------------------------------------------#
# make a TB model out of one of the solutions

cbr = candidatesv[1].apolarv[1]

# realize that if we only take intra-cell hoppings, the fitting will not converge
tbm = tb_hamiltonian(cbr, [[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1]]);

##-----------------------------------------------------------------------------------------#
# fit the TB model to the MPB results

# obtain the k-points and the spectrum
kp = irrfbz_path(sgnum, Rs)
kvs = interpolate(kp, 40)
ms = mpb.ModeSolver(;
    num_bands = 2,
    geometry_lattice = mp.Lattice(; basis1 = Rs[1], basis2 = Rs[2], basis3 = Rs[3]),
    geometry = pylist(geometry),
    k_points = pylist(map(k -> mp.Vector3(k...), kvs)),
)
ms.run()
freqs = pyconvert(Matrix{Float64}, ms.all_freqs)

μᴸ = tbm.N - size(freqs, 2) # number of longitudinal bands
cbrᴸ = candidatesv[1].longitudinal
μᴸ = occupation(cbrᴸ)
μᵀ = occupation(cbr) - μᴸ

ptbm_fit = photonic_fit(tbm, freqs[:, 1:μᵀ], kvs; verbose = true)
freqs_fit = spectrum(ptbm_fit, kvs; transform = energy2frequency)[:, μᴸ+1:end]

# ---------------------------------------------------------------------------------------- #
# plot the results

using GLMakie

plot(
    kvs,
    freqs,
    freqs_fit;
    color = [:blue, :red],
    linewidth = [3, 2],
    linestyle = [:solid, :dash],
)