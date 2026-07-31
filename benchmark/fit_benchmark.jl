# benchmark/fit_benchmark.jl
#
# Performance indicator for `fit` (SymmetricTightBindingOptimExt): fits synthetic reference
# spectra generated from models with known coefficients and reports, per scenario & draw,
# the fit time (ChairMarks `@b`: fastest sample; each sample re-seeds the RNG, so every
# sample performs identical work), allocations, the RMS energy error on the fitting
# k-points ("train") and on independent random k-points ("test"), and a success flag.
# All randomness is seeded, so results are reproducible and comparable across code changes:
# rerun after changing `fit`/gradients and diff the summary lines.
#
# Run from the repository root:
#
#     julia benchmark/fit_benchmark.jl

import Pkg
Pkg.activate(@__DIR__; io = devnull)

using Crystalline, SymmetricTightBinding, Optim, Chairmarks
using Random, Statistics, Printf, LinearAlgebra

# the near-singular-Hessian NaN workaround for `solve_tr_subproblem!` (hit by the SG 221
# scenario) now lives in `ext/optim_trsubproblem_patch.jl` and is applied automatically by
# `SymmetricTightBindingOptimExt` as soon as `Optim` is loaded above.

# include(joinpath(@__DIR__, "fit_mmd.jl")) # annealed-MMD prototype: superseded by
#                                           # basin-hopping `fit` (kept for reference)

# ---------------------------------------------------------------------------------------- #
# configuration

const N_DRAWS         = 5    # random true-coefficient draws per scenario
const N_K_RANDOM      = 100  # random fitting k-points (added to the high-symmetry ones);
                             # ~realistic sampling density (low values, ~10, overfit easily)
const N_K_TEST        = 50   # independent random k-points for the generalization error
const FIT_ATOL        = 1e-3 # `atol` handed to `fit`; also the noiseless success threshold
const BENCH_SECONDS   = 1.0  # `@b` time budget per draw: slow fits run once, fast fits
                             # report the fastest of several identical (re-seeded) runs

# fitters to compare: each is a function `(tbm, Em_r, ks; atol) -> ptbm` and manages its
# own trial/evaluation budget (each fitter's own default). `fit` (moment-seeded basin
# hopping) won the July 2026 global-search shootout vs. TikTak, θ-scheduled hopping,
# CMA-ES, and annealed-MMD prototypes; add experimental fitters here to compare against it
fitter_default(tbm, Em, ks; kws...) = fit(tbm, Em, ks; kws...)
const FITTERS = [
    "fit (basin-hop)" => fitter_default,
]

# ---------------------------------------------------------------------------------------- #
# scenarios
# NB: `@composite` is used at top-level only (it does not see local variables; cf. CLAUDE.md)

struct Scenario{D}
    name      :: String
    tbm       :: TightBindingModel{D} # the model to fit
    tbm_ref   :: TightBindingModel{D} # the model generating the reference (usually == tbm;
                                      # more terms = misspecified/truncated fit)
    cs_scales :: Vector{Float64}      # per-term std of the random reference coefficients
    hs_ks     :: Vector{Vector{Float64}} # high-symmetry k-points (stress degeneracies)
    σ_rel     :: Float64              # noise std relative to std of reference energies
end
function Scenario(name::String, tbm::TightBindingModel, hs_ks, σ_rel::Real)
    return Scenario(name, tbm, tbm, ones(length(tbm)), hs_ks, Float64(σ_rel))
end

println("building models (not part of the benchmark timings)…")
setup_t = @elapsed begin
    # graphene: (2b|A₁) EBR of plane group 17; 2 bands
    brs17  = calc_bandreps(17, Val(2))
    cbr17  = @composite brs17[5]
    tbm_nn  = tb_hamiltonian(cbr17, [[0, 0]])         # on-site + nearest neighbor
    tbm_nnn = tb_hamiltonian(cbr17, [[0, 0], [1, 0]]) # + next-nearest neighbor
    tbm_3rd = tb_hamiltonian(cbr17, [[0, 0], [1, 0], [1, 1]]) # + third-range hoppings
    hs17    = [[0.0, 0.0], [1/3, 1/3], [1/2, 0.0]]    # Γ, K (Dirac point), M

    # two-EBR graphene-lattice model: more bands, more terms, more crossings
    cbr17b  = @composite brs17[1] + brs17[5]
    tbm17b  = tb_hamiltonian(cbr17b, [[0, 0], [1, 0]])

    # SG 221: (3d|A₁g) ⊕ (3d|B₂g); 6 bands, 3D (the `fit` docstring example)
    brs221 = calc_bandreps(221)
    cbr221 = @composite brs221[1] + brs221[7]
    tbm221 = tb_hamiltonian(cbr221)
    hs221  = [[0.0, 0.0, 0.0], [0.0, 1/2, 0.0], [1/2, 1/2, 0.0], [1/2, 1/2, 1/2]] # Γ,X,M,R
end
@printf("   done (%.1f s)\n\n", setup_t)

# per-term coefficient scales for the longer-range graphene models: leading (on-site +
# shorter-range) terms O(1), longer-range terms small. NB: `tb_hamiltonian` does NOT order
# terms globally by hopping range (it is block-major, and only loosely range-sorted within a
# block); this trailing-terms-are-weak construction relies only on the weaker (verified)
# property that, for these single-EBR graphene models, `tbm_nnn`'s terms are a prefix of
# `tbm_3rd`'s and `tbm_4th`'s — so the trailing extra terms are the added longer-range ones
N_extra  = length(tbm_3rd) - length(tbm_nnn)
N_extra4 = length(tbm_4th) - length(tbm_nnn)
scales_hier  = vcat(ones(length(tbm_nnn)), fill(0.2, N_extra)) # realistic hopping decay
scales_trunc = vcat(ones(length(tbm_nnn)), fill(0.1, N_extra)) # small omitted terms
scales_trunc4 = vcat(ones(length(tbm_nnn)), [0.2^r for r in 1:N_extra4]) # strongly decaying tail

scenarios = [
    Scenario("graphene NN",            tbm_nn,  hs17,  0.0),
    Scenario("graphene NNN",           tbm_nnn, hs17,  0.0),
    Scenario("graphene NNN, 2% noise", tbm_nnn, hs17,  0.02),
    Scenario("SG 221 (6 bands, 3D)",   tbm221,  hs221, 0.0),
    # exact-model fit with realistically decaying (hierarchical) coefficient scales
    Scenario("graphene 3-range, hierarchical", tbm_3rd, tbm_3rd, scales_hier, hs17, 0.0),
    # misspecified fit: reference has N terms, fit model only M < N (omitted terms small);
    # no exact solution exists — success means beating the naive-truncation error Δ
    Scenario("graphene truncated (M<N)", tbm_nnn, tbm_3rd, scales_trunc, hs17, 0.0),
    # larger misspecified fit with a strongly range-decaying reference (4 hopping ranges, fit
    # keeps only the 4 short-range terms): the frustrated, hierarchical regime that hopping-
    # range continuation targets — cf. `fit_continuation`
    Scenario("graphene truncated 4-range (M<N)", tbm_nnn, tbm_4th, scales_trunc4, hs17, 0.0),
    # two EBRs: more bands, more terms, more crossings — stresses the multi-start search
    Scenario("pg 17, two EBRs",        tbm17b,  hs17,  0.0),
]

# ---------------------------------------------------------------------------------------- #
# benchmark driver

rms(A) = sqrt(mean(abs2, A))

function run_scenario(s::Scenario, scenario_idx::Integer, fitter)
    D  = length(first(s.hs_ks))
    Nᶜ = length(s.tbm)
    Nᵇ = s.tbm.N
    Nᵏ = length(s.hs_ks) + N_K_RANDOM
    truncated = length(s.tbm_ref) > Nᶜ # misspecified fit: reference has more terms
    @printf("== %s  (%d terms%s, %d bands, D=%d, %d k-points%s) ==\n",
            s.name, Nᶜ, truncated ? " of $(length(s.tbm_ref))" : "", Nᵇ, D, Nᵏ,
            s.σ_rel > 0 ? ", noisy" : "")
    println("   draw   time (s)   alloc (MiB)   RMS train   RMS test    success")

    results = NamedTuple[]
    for draw in 1:N_DRAWS
        rng = Xoshiro(1000scenario_idx + draw) # data seed: fixed per (scenario, draw)
        cs_true   = s.cs_scales .* randn(rng, length(s.tbm_ref))
        ptbm_true = s.tbm_ref(cs_true)
        ks   = vcat(s.hs_ks, [rand(rng, D) for _ in 1:N_K_RANDOM])
        Em_r = spectrum(ptbm_true, ks)
        σ = s.σ_rel * std(Em_r)
        σ > 0 && (Em_r = Em_r .+ σ .* randn(rng, size(Em_r)))

        # naive-truncation model error Δ: reference bands vs. reference with the omitted
        # (trailing; cf. `scales_trunc` ordering note) terms zeroed. A misspecified fit
        # cannot be exact, but it can absorb part of the omitted terms — so it should at
        # least beat Δ
        Δ = if truncated
            cs_naive = copy(cs_true)
            cs_naive[(Nᶜ+1):end] .= 0.0
            rms(Em_r .- spectrum(s.tbm_ref(cs_naive), ks))
        else
            0.0
        end

        # benchmark the fit; the RNG re-seeding makes every sample perform identical work,
        # so the fastest sample is a clean timing of a fully deterministic computation
        fitted = Ref{Any}()
        bench = @b begin
            Random.seed!(1234 + draw)
            fitted[] = fitter(s.tbm, Em_r, ks;
                              atol = max(FIT_ATOL, σ, Δ/2)) # don't demand the unattainable
        end seconds = BENCH_SECONDS
        ptbm_fit = fitted[]
        t   = bench.time
        mib = bench.bytes / 2^20

        rms_train = rms(spectrum(ptbm_fit, ks) .- Em_r)
        ks_test   = [rand(rng, D) for _ in 1:N_K_TEST]
        rms_test  = rms(spectrum(ptbm_fit, ks_test) .- spectrum(ptbm_true, ks_test))
        # noiseless & exact: train RMS ≤ FIT_ATOL; noisy: test RMS within 2× noise floor;
        # truncated: at least as good as naive truncation
        success = σ > 0    ? (rms_test ≤ 2σ) :
                  truncated ? (rms_train ≤ Δ) : (rms_train ≤ FIT_ATOL)

        note = truncated ? (@sprintf "   Δ_naive=%.1e" Δ) : ""
        @printf("   %4d   %8.3f   %11.1f   %9.2e   %9.2e   %s%s\n",
                draw, t, mib, rms_train, rms_test, success ? "✓" : "✗", note)
        push!(results, (; draw, t, mib, rms_train, rms_test, success))
    end

    n_success = count(r -> r.success, results)
    med_t     = median(r.t for r in results)
    med_mib   = median(r.mib for r in results)
    med_test  = median(r.rms_test for r in results)
    @printf("   summary: %d/%d success | median time %.3f s | median alloc %.1f MiB | median test-RMS %.2e\n\n",
            n_success, N_DRAWS, med_t, med_mib, med_test)
    return results
end

# ---------------------------------------------------------------------------------------- #
# warm-up (compile the 2D & 3D fit paths so first-draw timings are not dominated by JIT)

print("warming up (compilation)… ")
warmup_t = @elapsed for (tbm, D) in ((tbm_nn, 2), (tbm221, 3))
    Random.seed!(1)
    ks = [rand(D) for _ in 1:3]
    Em = spectrum(tbm(randn(length(tbm))), ks)
    for (_, fitter) in FITTERS
        fitter(tbm, Em, ks; atol = 1e6) # huge atol → early return after the first trial
    end
end
@printf("done (%.1f s)\n\n", warmup_t)

# ---------------------------------------------------------------------------------------- #
# run & report

gitrev = try
    root  = dirname(@__DIR__)
    dirty = !isempty(read(`git -C $root status --porcelain`, String))
    readchomp(`git -C $root rev-parse --short HEAD`) * (dirty ? "-dirty" : "")
catch
    "unknown"
end
println("fit benchmark @ SymmetricTightBinding ", gitrev, " | julia ", VERSION, "\n")

all_results = map(FITTERS) do (name, fitter)
    printstyled("―― fitter: ", name, " ――\n\n"; bold = true)
    fitter_results = [run_scenario(s, i, fitter) for (i, s) in enumerate(scenarios)]
    n_success = sum(count(r -> r.success, rs) for rs in fitter_results)
    n_total   = sum(length, fitter_results)
    t_total   = sum(sum(r.t for r in rs) for rs in fitter_results)
    @printf("―― TOTAL (%s): %d/%d success | total fit time %.1f s ――\n\n",
            name, n_success, n_total, t_total)
    fitter_results
end

# comparison table: rows = fitters, columns = scenarios; entries = success & median time
println("== COMPARISON (success | median fit time) ==")
@printf("%-18s", "")
foreach(s -> @printf(" | %-22s", s.name), scenarios)
println()
for ((name, _), fitter_results) in zip(FITTERS, all_results)
    @printf("%-18s", name)
    for rs in fitter_results
        @printf(" | %d/%d %8.3f s%7s", count(r -> r.success, rs), N_DRAWS,
                median(r.t for r in rs), "")
    end
    println()
end
