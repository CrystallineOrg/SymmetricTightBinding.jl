using Pkg
Pkg.activate(@__DIR__)

### necessary packages
using TETB
using Crystalline
using TETB.PythonCall: pylist
using TETB.SymmetricTightBinding: ReciprocalPointLike
using Optim

### construct the structure under study

R1 = 0.2 #cylinder radius
mat = mp.Medium(; epsilon = 12)
geometry = [
    mp.Cylinder(;
        radius = R1,
        center = [0, 0, 0],
        axis = [0, 0, 1],
        height = 1,
        material = mat,
    ),
    mp.Cylinder(;
        radius = R1,
        center = [0, 0, 0],
        axis = [0, 1, 0],
        height = 1,
        material = mat,
    ),
    mp.Cylinder(;
        radius = R1,
        center = [0, 0, 0],
        axis = [1, 0, 0],
        height = 1,
        material = mat,
    ),
]

### solve the system
ms = mpb.ModeSolver(;
    num_bands = 8,
    geometry_lattice = mp.Lattice(;
        basis1 = [1, 0, 0],
        basis2 = [0, 1, 0],
        basis3 = [0, 0, 1],
        size = [1, 1, 1],
    ),
    geometry = pylist(geometry),
    resolution = 16,
)
ms.init_params(; p = mp.ALL, reset_fields = true)

### obtain the symmetry vectors of the bands computed above
sgnum = 221
symvecs = obtain_symmetry_vectors(ms, sgnum)

m = symvecs[1] # pick the 2 lower bands

### obtain an EBR decomposition with at least one additional band
μᴸ = 1
brs = calc_bandreps(sgnum)
candidatesv = find_bandrep_decompositions(m, brs; μᴸ_min = μᴸ)

##-----------------------------------------------------------------------------------------#
# make a TB model out of one of the solutions

cbr = candidatesv[1].apolarv[1]

tbm = tb_hamiltonian(cbr, [[0, 0, 0]])

##-----------------------------------------------------------------------------------------#
# fit the TB model to the MPB results

# obtain the k-points and the spectrum
# TODO: ↓ maybe we can perform this step using Brillouin instead of MPB interpolation
k_points = [
    mp.Vector3(),  # Gamma
    mp.Vector3(0, 0.5, 0),  # X
    mp.Vector3(0.5, 0.5, 0),  # M
    mp.Vector3(),  # Gamma
    mp.Vector3(0, 0.5, 0.5),  # R
    mp.Vector3(0, 0.5, 0),  # X
]
k_interp = 5 # number of k-points
k_points = mp.interpolate(k_interp, k_points)
ms = mpb.ModeSolver(;
    num_bands = 2,
    geometry_lattice = mp.Lattice(;
        basis1 = [1, 0, 0],
        basis2 = [0, 1, 0],
        basis3 = [0, 0, 1],
        size = [1, 1, 1],
    ),
    geometry = pylist(geometry),
    k_points = k_points,
)
ms.run()
freqs = ms.all_freqs

# transform the interesting variables to Julia types
Em_r = TETB.PythonCall.pyconvert(Matrix, freqs)
ks = TETB.PythonCall.pyconvert(Vector{ReciprocalPoint{3}}, k_points)

# ptbm_fit = fit(tbm, Em_r, ks) # TODO: it does not converge...
# Em_fitted = spectrum(ptbm_fit, ks)

# ---------------------------------------------------------------------------------------- #
# plot the results

using GLMakie

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "k", ylabel = "Frequency (c/2πa)")

for i in axes(Em_r, 2)
    lines!(ax, Em_r[:, i]; linestyle = :dash, color = :black, linewidth = 1)
end

fig