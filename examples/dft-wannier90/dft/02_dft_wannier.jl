# Step 2 (DFT environment, `examples/dft-wannier90/dft/Project.toml`)
#
# Read `out/<case>/setup.toml` (from step 1), then:
#   1. run a self-consistent DFT calculation (DFTK.jl),
#   2. compute the DFT band structure along a high-symmetry path → `bands_dft.{dat,toml}`,
#   3. Wannierize (Wannier90, via DFTK's interface and `wannier90_jll`) with the initial
#      projections of `setup.toml` → `wannier/wann_hr.dat` (and `_wsvec.dat`, `.wout`, ...).
# All energies written are relative to an energy zero: the valence-band maximum (VBM) or the
# Fermi level (`setup.toml`: `dft.reference`); `_hr.dat` is in DFTK's absolute energies.
#
# Usage: julia --project=dft dft/02_dft_wannier.jl [case] [mode]
#   `mode = "bands"`: steps 1–2 only (e.g., to inspect the bands before choosing windows)

using DFTK, PseudoPotentialData, wannier90_jll
using Unitful, UnitfulAtomic
using TOML, Printf, LinearAlgebra

case = get(ARGS, 1, "si-valence")
mode = get(ARGS, 2, "all")
outdir = joinpath(@__DIR__, "..", "out", case)
setup = TOML.parsefile(joinpath(outdir, "setup.toml"))
dft, wann = setup["dft"], setup["wannier"]
eV(ε) = ustrip(auconvert(u"eV", ε)) # Hartree → eV

# ---------------------------------------------------------------------------------------- #
# 1. SCF

lattice = austrip.(reduce(hcat, setup["lattice_angstrom"]) .* u"Å") # columns, in bohr
family = PseudoFamily(dft["pseudo_family"])
# NB: identical atoms must share one element object; otherwise, DFTK (spglib) treats them
#     as distinct species and finds too few symmetries (e.g., Pmmm instead of Pm-3m)
elements = Dict(s => ElementPsp(Symbol(s), family) for s in unique(setup["atoms"]["species"]))
atoms = [elements[s] for s in setup["atoms"]["species"]]
positions = [Vector{Float64}(r) for r in setup["atoms"]["positions"]]
smearing_kws = haskey(dft, "temperature") ?
    (; temperature = dft["temperature"], smearing = DFTK.Smearing.Gaussian()) : (;)
model = model_DFT(lattice, atoms, positions; functionals = LDA(), smearing_kws...)

n_bands_w90 = haskey(wann, "bands") ? wann["bands"][2] + 1 : wann["n_bands"]
n_bands_conv = max(n_bands_w90, dft["n_bands"])
basis = PlaneWaveBasis(model; Ecut = dft["Ecut"], kgrid = dft["kgrid"])
nbandsalg = AdaptiveBands(model; n_bands_converge = n_bands_conv)
scfres = self_consistent_field(basis; nbandsalg, tol = 1e-8)

# energy zero
εs = [eV.(εk) for εk in scfres.eigenvalues] # per irreducible k-point, eV
εF = eV(scfres.εF)
E₀ = if dft["reference"] == "vbm" # insulators: top of the occupied bands
    n_occ = count(>(1e-3), scfres.occupation[1])
    maximum(εk[n_occ] for εk in εs)
elseif dft["reference"] == "fermi"
    εF
else
    error("unknown energy reference $(dft["reference"])")
end
@printf "energy zero (%s) = %.4f eV; Fermi level = %.4f eV\n" dft["reference"] E₀ εF
# band ranges, relative to the energy zero
for n in 1:n_bands_conv
    @printf "  band %2d: [% 8.3f, % 8.3f] eV\n" n minimum(εk[n] for εk in εs)-E₀ maximum(εk[n] for εk in εs)-E₀
end

# ---------------------------------------------------------------------------------------- #
# 2. DFT band structure along a high-symmetry path (DFTK's `irrfbz_path`, via Brillouin.jl);
#    k-points are in reduced coordinates of the reciprocal basis of `lattice`, i.e.,
#    Crystalline's primitive reciprocal basis. For 2D cases, the path lies in the k_z = 0
#    plane (via the plane group number)

kpath = setup["dim"] == 2 ?
    irrfbz_path(model; dim = 2, space_group_number = setup["sgnum"]) : irrfbz_path(model)
bands = compute_bands(scfres, kpath; n_bands = dft["n_bands"],
                      kline_density = dft["kline_density"], tol = 1e-8)
kinter = bands.kinter
open(joinpath(outdir, "bands_dft.dat"), "w") do io
    println(io, "# branch  k₁ k₂ k₃ (reduced, primitive)  E₁ … E_n (eV, relative to the energy zero)")
    ik = 0
    for (b, kbranch) in enumerate(kinter.kpaths), _ in kbranch
        ik += 1
        k = bands.basis.kpoints[ik].coordinate
        @printf io "%d  % .10f % .10f % .10f " b k...
        for ε in eV.(bands.eigenvalues[ik])
            @printf io " % .8f" ε - E₀
        end
        println(io)
    end
end
# high-symmetry labels: (branch, index within branch, label)
labels = sort!([Dict("branch" => b, "index" => i, "label" => string(l))
                for (b, d) in enumerate(kinter.labels) for (i, l) in d];
               by = d -> (d["branch"], d["index"]))
open(io -> TOML.print(io, Dict("reference" => dft["reference"], "reference_eV" => E₀,
                               "fermi_eV" => εF, "labels" => labels)),
     joinpath(outdir, "bands_dft.toml"), "w")
mode == "bands" && exit()

# ---------------------------------------------------------------------------------------- #
# 3. Wannierization

projections = map(setup["projections"]) do p
    # hydrogenic, nodeless (Wannier90's `r = 1`); α = 1/Å (Wannier90's default `zona`)
    DFTK.HydrogenicWannierProjection(Vector{Float64}(p["center"]), 1, p["l"], p["m"],
                                     1 / austrip(1u"Å"))
end

w90kws = Dict{Symbol, Any}(
    :num_iter => wann["num_iter"],
    :write_hr => true,         # → `_hr.dat`
    :write_xyz => true,        # → `_centres.xyz`
    :translate_home_cell => false,
    :use_ws_distance => true,  # minimal-image hoppings (Wannier90 default; → `_wsvec.dat`)
    :bands_plot => true,       # Wannier90's own interpolation, as a cross-check
)
if haskey(wann, "bands") # isolated bands b₁…b₂: windows in the gaps around them
    b₁, b₂ = wann["bands"]
    b₂ - b₁ + 1 == wann["n_wannier"] || error("number of bands ≠ number of projections")
    if b₁ > 1
        εmax(n) = maximum(εk[n] for εk in εs)
        εmin(n) = minimum(εk[n] for εk in εs)
        (εmax(b₁ - 1) < εmin(b₁) && εmax(b₂) < εmin(b₂ + 1)) ||
            error("bands $b₁:$b₂ are not isolated")
        lo, hi = (εmax(b₁ - 1) + εmin(b₁)) / 2, (εmax(b₂) + εmin(b₂ + 1)) / 2
        @printf "windows (isolated bands %d:%d): [%.3f, %.3f] eV\n" b₁ b₂ lo-E₀ hi-E₀
        merge!(w90kws, Dict(:dis_win_min => lo, :dis_win_max => hi,
                            :dis_froz_min => lo, :dis_froz_max => hi))
    else
        global n_bands_w90 = b₂ # lowest bands: no windows needed
    end
else # entangled bands: windows relative to the energy zero
    for key in ("dis_win_min", "dis_win_max", "dis_froz_min", "dis_froz_max")
        haskey(wann, key) && (w90kws[Symbol(key)] = E₀ + wann[key])
    end
    w90kws[:dis_num_iter] = wann["dis_num_iter"]
end

fileprefix = joinpath(outdir, "wannier", "wann")
run_wannier90(scfres; fileprefix, n_bands = n_bands_w90, n_wannier = wann["n_wannier"],
              projections, w90kws...)

println("done: ", relpath(fileprefix * "_hr.dat", joinpath(@__DIR__, "..")))
