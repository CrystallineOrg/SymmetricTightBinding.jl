# Step 3 (analysis environment, `examples/dft-wannier90/Project.toml`)
#
# Read the Wannier90 `_hr.dat` file of step 2, build the symmetric tight-binding model `tbm`
# of the band representation of step 1 over all translations of the Wannier Hamiltonian,
# and fit its coefficients to the Wannier90 hopping matrices (`fit_wannier90`; exact if the
# Wannier Hamiltonian is symmetric). Then extract a *compact* model: rank the hopping orbits
# by the norm of their amplitudes and keep the fewest, strongest orbits whose bands deviate
# from the Wannier90 bands by at most `tol` (RMS, along the DFT k-path). Finally, compare
# band structures: DFT vs. Wannier90 vs. the compact model. Writes `out/<case>/bands.png`,
# `out/<case>/fit.toml`.
#
# Usage: julia --project=. 03_fit.jl [case] [tol]        (`tol` in eV; default 0.05)

using Crystalline, SymmetricTightBinding
using Crystalline: constant
using SymmetricTightBinding: read_wannier90, fit_wannier90, _orbit_term_ranges,
                             _fit_hopping_orbit
using TOML, DelimitedFiles, LinearAlgebra, Printf
using CairoMakie

include("common.jl")

case = get(ARGS, 1, "si-valence")
tol = parse(Float64, get(ARGS, 2, "0.05"))
outdir = joinpath(@__DIR__, "out", case)
setup = TOML.parsefile(joinpath(outdir, "setup.toml"))
bands_info = TOML.parsefile(joinpath(outdir, "bands_dft.toml"))
spec = CASES[case]
D = spec.dim
cbr = composite_from_labels(bandreps(spec.sgnum, D), spec.brlabels)
positions = orbital_positions(cbr)
lattice = reduce(hcat, setup["lattice_angstrom"])[1:D, 1:D] # columns, Å
E₀ = bands_info["reference_eV"] # energy zero (VBM or Fermi level)
rms(x) = sqrt(sum(abs2, x) / length(x))

# ---------------------------------------------------------------------------------------- #
# 1. read the Wannier90 model; energies relative to the energy zero; basis alignment

prefix = joinpath(outdir, "wannier", "wann")
w90tbm = read_wannier90(prefix, positions) # `wann_hr.dat` & `wann_wsvec.dat`
for (R, H) in zip(w90tbm.translations, w90tbm.hopping_matrices)
    iszero(R) && (H[diagind(H)] .-= E₀) # shift on-site energies
end
println("Wannier90 model: ", orbital_count(w90tbm), " orbitals, ",
        length(w90tbm.translations), " translations")

# cross-check against Wannier90's own band interpolation (`_band.kpt` & `_band.dat`)
let kpt = readdlm(prefix * "_band.kpt"; skipstart = 1)
    ks_w90 = [kpt[i, 1:D] for i in axes(kpt, 1)]
    E_w90 = reshape(readdlm(prefix * "_band.dat")[:, 2], length(ks_w90), :) .- E₀
    E_ours = stack(k -> eigvals(Hermitian(w90tbm(k))), ks_w90; dims = 1)
    @printf "cross-check vs. Wannier90's `_band.dat`: max |ΔE| = %.2e eV\n" maximum(abs, E_ours - E_w90)
end

# express the hopping matrices in the local orbital basis of the tight-binding model (see
# `basis_alignment` in `common.jl`)
lms = [(p["l"], p["m"]) for p in setup["projections"]]
U = basis_alignment(cbr, lms, lattice)
if !(U ≈ I)
    println("basis alignment |stbₐ⟩ = ∑_b |w90_b⟩ U[b,a], with U = ")
    show(stdout, "text/plain", round.(U; digits = 4)); println()
end
w90tbm = align_basis(w90tbm, U)

# ---------------------------------------------------------------------------------------- #
# 2. the full symmetric model, over all translations of the Wannier Hamiltonian, & its fit

t_tbm = @elapsed tbm = tb_hamiltonian(cbr, [Vector(R) for R in w90tbm.translations])
# the per-orbit fits (as in `fit_wannier90`), with their residuals: for diagnostics, and
# to rank the orbits for the compact model below
w90_hoppings = Dict(zip(w90tbm.translations, w90tbm.hopping_matrices))
orbit_fits = [(; terms, _fit_hopping_orbit(tbm, terms, w90_hoppings)...)
              for terms in _orbit_term_ranges(tbm)]
ptbm = fit_wannier90(tbm, w90tbm)
@printf "full model: %d terms in %d hopping orbits (tb_hamiltonian: %.1f s)\n" length(tbm) length(orbit_fits) t_tbm
@printf "fit: max residual over orbits %.1e eV; total relative residual %.1e\n" maximum(o -> o.residual, orbit_fits) norm(o.residual for o in orbit_fits) / norm(o.w90_norm for o in orbit_fits)

# ---------------------------------------------------------------------------------------- #
# 3. band structures along the DFT k-path

dat = readdlm(joinpath(outdir, "bands_dft.dat"); comments = true)
branches = Int.(dat[:, 1])
ks = [dat[i, 2:1+D] for i in axes(dat, 1)]
E_dft = dat[:, 5:end]
E_w90 = stack(k -> eigvals(Hermitian(w90tbm(k))), ks; dims = 1)
nw = orbital_count(w90tbm)
@printf "full model vs. Wannier90: max |ΔE| = %.1e eV\n" maximum(abs, spectrum(ptbm, ks) - E_w90)

# error vs. DFT, over the DFT states the Wannier functions must reproduce: an isolated set of
# bands (by band index); or, for entangled bands, the states in the frozen window (compared
# with the nearest Wannier/model band)
wann = setup["wannier"]
function Δdft(E)
    if haskey(wann, "bands")
        b₁, b₂ = wann["bands"]
        return vec(E - E_dft[:, b₁:b₂])
    else
        lo, hi = get(wann, "dis_froz_min", -Inf), wann["dis_froz_max"]
        return [minimum(e -> abs(e - ε), E[i, :]) for i in axes(E_dft, 1)
                for ε in E_dft[i, :] if lo ≤ ε ≤ hi]
    end
end

# ---------------------------------------------------------------------------------------- #
# 4. compact model: the strongest hopping orbits, ranked by the norm of their amplitudes

orbit_δ(o) = constant(first(orbit(tbm.terms[first(o.terms)].block.h_orbit)))
orbit_block(o) = tbm.terms[first(o.terms)].block_ij
ranked = sort(orbit_fits; by = o -> -o.w90_norm)
function compact_model(n)
    keep = sort!(reduce(vcat, collect(o.terms) for o in ranked[1:n]))
    return tbm[keep](ptbm.cs[keep]), keep
end
println("\ncompact models (strongest n orbits):")
println("  orbits  terms  max |δ| (Å)   vs. Wannier90: RMS / max (eV)   vs. DFT: RMS / max (eV)")
n_c = nothing
for n in 1:length(ranked)
    pc, _ = compact_model(n)
    E = spectrum(pc, ks)
    err_w90, err_dft = E - E_w90, Δdft(E)
    if n ≤ 12 || rms(err_w90) ≤ tol || n == length(ranked) || n % 10 == 0
        @printf "  %5d  %5d   %8.3f        %.3f / %.3f                  %.3f / %.3f\n" n length(pc.cs) maximum(o -> norm(lattice * orbit_δ(o)), ranked[1:n]) rms(err_w90) maximum(abs, err_w90) rms(err_dft) maximum(abs, err_dft)
    end
    if rms(err_w90) ≤ tol
        global n_c = n
        break
    end
end
ptbm_c, keep = compact_model(n_c)
E_c = spectrum(ptbm_c, ks)
@printf "\ncompact model: %d orbits, %d terms (RMS ≤ %.3f eV vs. Wannier90)\n" n_c length(keep) tol
println("  block    δ (reduced)                     |δ| (Å)   cs")
for o in sort(ranked[1:n_c]; by = o -> (norm(lattice * orbit_δ(o)), orbit_block(o)))
    @printf "  %-8s %-30s %7.3f   %s\n" orbit_block(o) round.(orbit_δ(o); digits = 4) norm(lattice * orbit_δ(o)) join((@sprintf("% .4f", c) for c in o.c), " ")
end
println()
show(stdout, MIME"text/plain"(), ptbm_c); println()

# ---------------------------------------------------------------------------------------- #
# 5. plot: DFT vs. Wannier90 vs. the compact model

B = 2π * inv(lattice)' # reciprocal basis (columns), 1/Å
xs = zeros(length(ks)) # x-axis: cumulative Cartesian distance (1/Å); branches concatenated
for i in 2:length(ks)
    xs[i] = xs[i-1] + (branches[i] == branches[i-1] ? norm(B * (ks[i] - ks[i-1])) : 0.0)
end
branch_starts = [findfirst(==(b), branches) for b in unique(branches)]
ticks = Dict{Float64, String}()
for l in bands_info["labels"]
    x = xs[branch_starts[l["branch"]] + l["index"] - 1]
    ticks[x] = haskey(ticks, x) && ticks[x] != l["label"] ? ticks[x] * "|" * l["label"] : l["label"]
end
tick_xs = sort!(collect(keys(ticks)))

E_lims = (minimum(E_w90) - 1.5, maximum(E_w90) + 1.5)
fig = Figure(size = (720, 520))
ax = Axis(fig[1, 1]; ylabel = "Energy − $(bands_info["reference"] == "vbm" ? "VBM" : "E_F") (eV)",
          xticks = (tick_xs, [ticks[x] for x in tick_xs]),
          title = "$(spec.description)\n$(cbr): DFT vs. Wannier90 vs. compact model",
          titlesize = 13)
vlines!(ax, tick_xs; color = :gray80, linewidth = 0.8)
hlines!(ax, [0.0]; color = :gray80, linewidth = 0.8, linestyle = :dash)
for b in unique(branches) # avoid connecting lines across branch jumps
    idx = findall(==(b), branches)
    for n in axes(E_dft, 2)
        lines!(ax, xs[idx], E_dft[idx, n]; color = :black, linewidth = 2.5,
               label = (b, n) == (1, 1) ? "DFT" : nothing)
    end
    for n in 1:nw
        lines!(ax, xs[idx], E_w90[idx, n]; color = :dodgerblue, linewidth = 1.5,
               label = (b, n) == (1, 1) ? "Wannier90" : nothing)
        lines!(ax, xs[idx], E_c[idx, n]; color = :orangered, linewidth = 1.5,
               linestyle = :dash,
               label = (b, n) == (1, 1) ? "compact model ($(length(keep)) terms)" : nothing)
    end
end
xlims!(ax, first(xs), last(xs))
ylims!(ax, E_lims...)
Legend(fig[2, 1], ax; orientation = :horizontal, framevisible = false)
save(joinpath(outdir, "bands.png"), fig; px_per_unit = 2)
println("\nwrote ", relpath(joinpath(outdir, "bands.png"), @__DIR__))

open(joinpath(outdir, "fit.toml"), "w") do io
    TOML.print(io, Dict(
        "cbr" => string(cbr), "tol_eV" => tol,
        "full" => Dict("n_terms" => length(tbm), "n_orbits" => length(orbit_fits),
                       "max_orbit_residual_eV" => maximum(o -> o.residual, orbit_fits)),
        "compact" => Dict("n_orbits" => n_c, "terms" => keep, "cs" => ptbm_c.cs,
                          "rms_vs_w90_eV" => rms(E_c - E_w90),
                          "rms_vs_dft_eV" => rms(Δdft(E_c))),
        "rms_w90_vs_dft_eV" => rms(Δdft(E_w90))))
end
