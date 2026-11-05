# Step 1 (analysis environment, `examples/dft-wannier90/Project.toml`)
#
# Fix the crystal structure and the band representation (BR) of the bands to Wannierize,
# and write everything the DFT + Wannier90 step needs to `out/<case>/setup.toml`: the
# lattice and atoms (in Crystalline's primitive setting, so that the DFT, Wannier90, and
# SymmetricTightBinding all refer to the same basis), and one Wannier90 initial projection
# per tight-binding orbital, centred at the orbital's Wyckoff position and listed in
# SymmetricTightBinding's orbital order (so that the Wannier functions come out in the
# same order as the orbitals of the tight-binding model).
#
# Two-dimensional cases are embedded in 3D for the DFT: the third lattice vector is
# perpendicular to the plane, of length `vacuum`, with all atoms in the plane z = 0.
#
# Usage: julia --project=. 01_setup.jl [case]          (cases: see `CASES` in `common.jl`)

using Crystalline, SymmetricTightBinding
using TOML

include("common.jl")

case = get(ARGS, 1, "si-valence")
spec = CASES[case]
D = spec.dim
brs = bandreps(spec.sgnum, D)
cbr = composite_from_labels(brs, spec.brlabels)
pad3(v) = [v; zeros(3 - length(v))]

# primitive lattice (Crystalline's setting), in Å; columns are the lattice vectors
Rs = directbasis(spec.sgnum, Val(D)) # conventional; a = 1 for cubic & hexagonal systems
pRs = primitivize(Rs, centering(spec.sgnum, D))
lattice = [pad3(spec.a * pRs[i]) for i in 1:D]
D == 2 && push!(lattice, [0.0, 0.0, spec.vacuum])
D == 3 || D == 2 || error("D = $D not set up")

# atoms, at the sites of their Wyckoff positions
species = String[]
atom_positions = Vector{Float64}[]
for (el, wp_label) in spec.atoms, r in wyckoff_sites(brs, wp_label)
    push!(species, el)
    push!(atom_positions, pad3(collect(r)))
end

# projections: one per orbital, in SymmetricTightBinding's orbital order
positions = orbital_positions(cbr)
projections = [Dict("center" => pad3(collect(r)), "l" => l, "m" => m)
               for (r, (l, m)) in zip(positions, projection_angular_parts(cbr, spec.angular))]

setup = Dict(
    "case" => case,
    "description" => spec.description,
    "sgnum" => spec.sgnum,
    "dim" => D,
    "cbr" => string(cbr),
    "lattice_angstrom" => lattice,
    "atoms" => Dict("species" => species, "positions" => atom_positions),
    "dft" => spec.dft,
    "wannier" => merge(spec.wannier, Dict("n_wannier" => length(positions))),
    "projections" => projections,
)

outdir = mkpath(joinpath(@__DIR__, "out", case))
open(joinpath(outdir, "setup.toml"), "w") do io
    TOML.print(io, setup; sorted = true)
end
println("wrote ", relpath(joinpath(outdir, "setup.toml"), @__DIR__), " for ", cbr, " (",
        length(positions), " orbitals)")
