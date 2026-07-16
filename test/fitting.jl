using Test
using SymmetricTightBinding
using SymmetricTightBinding: spectralmoments
using Crystalline
using Optim
using Random
using LinearAlgebra: Hermitian, eigen, tr

# Tests for the Optim.jl fitting extension (`fit`, `multistart_fit`, `make_objective`,
# `spectralmoments`) and the `TightBindingCache` it is built on. Kept deliberately fast:
# small 2D graphene models plus one small 3D model, few k-points, and exactly-representable
# reference spectra (so the multi-start search early-returns quickly).

rms(A) = sqrt(sum(abs2, A) / length(A))

@testset "Fitting (Optim extension)" begin
    # --- shared small models -------------------------------------------------------------
    brs = calc_bandreps(17, Val(2))
    cbr = @composite brs[5]                                 # graphene honeycomb, (2b|A₁)
    tbm_nn  = tb_hamiltonian(cbr, [[0, 0], [1, 0]])         # on-site + nearest neighbour
    tbm_nnn = tb_hamiltonian(cbr, [[0, 0], [1, 0], [1, 1]]) # + next-nearest neighbour

    # a 3D two-EBR model (SG 221, (3d|A₁g) ⊕ (3d|B₂g)): 6 bands, 4 terms. Terms 2 & 4 are the
    # inter-EBR coupling blocks and are purely off-diagonal (traceless at every k) — used by
    # both the 3D recovery test and the rank-deficient moment-seed regression below.
    brs221 = calc_bandreps(221)
    cbr221 = @composite brs221[1] + brs221[7]
    tbm221 = tb_hamiltonian(cbr221)

    # a small, fixed set of k-points: three high-symmetry points plus random filler. Enough
    # to over-determine the (few) hopping amplitudes without making the tests slow.
    ks = let
        Random.seed!(0)
        vcat([[0.0, 0.0], [1/3, 1/3], [1/2, 0.0]], [rand(2) for _ in 1:12])
    end

    # -------------------------------------------------------------------------------------
    @testset "TightBindingCache matches direct evaluation" begin
        Random.seed!(1)
        cs = randn(length(tbm_nnn))
        ptbm = tbm_nnn(cs)
        cache = TightBindingCache(tbm_nnn, ks)

        # assembly parity: cache(cs, κ) == ∑ᵢ csᵢ hᵢ(ks[κ]), and eigenvalues match `spectrum`
        Em = spectrum(ptbm, ks)
        for κ in eachindex(ks)
            H_cache = Matrix(cache(cs, κ)) # peel the Hermitian wrapper
            H_sum = sum(c .* tbm_nnn[i](ks[κ]) for (i, c) in enumerate(cs))
            @test H_cache ≈ H_sum
            @test eigen(Hermitian(copy(H_cache))).values ≈ Em[κ, :]
        end

        # gradient parity: cached vs uncached Feynman–Hellmann coefficient gradient, using
        # the same eigensolution for both (no Bloch phase, as the fitting loss uses)
        for κ in eachindex(ks)
            Es, us = eigen(Hermitian(Matrix(cache(cs, κ))))
            ∇c = collect(collect(g) for g in energy_gradient_wrt_hopping(cache, κ, (Es, us)))
            ∇d = collect(collect(g) for g in energy_gradient_wrt_hopping(ptbm, ks[κ], (Es, us)))
            @test ∇c ≈ ∇d
        end
    end

    # -------------------------------------------------------------------------------------
    @testset "make_objective / spectralmoments routes agree" begin
        Random.seed!(2)
        cs_r = randn(length(tbm_nn))
        Em_r = spectrum(tbm_nn(cs_r), ks)
        cache = TightBindingCache(tbm_nn, ks)

        # the cache-based and (tbm, ks)-based constructors must build the same objects
        m_cache = spectralmoments(cache, Em_r)
        m_tbm   = spectralmoments(tbm_nn, Em_r, ks)
        @test m_cache.c₀ ≈ m_tbm.c₀
        @test m_cache.ρ ≈ m_tbm.ρ
        @test m_cache.ρ > 0 && all(>(0), m_cache.ρs)

        # both `make_objective` routes are usable and produce a fit of the same quality
        for obj in (make_objective(cache, Em_r), make_objective(tbm_nn, Em_r, ks))
            best_cs, best_loss = multistart_fit(obj, m_cache; max_multistarts = 40)
            @test best_loss < 1e-8
            @test rms(spectrum(tbm_nn(best_cs), ks) - Em_r) < 1e-5
        end
    end

    # -------------------------------------------------------------------------------------
    @testset "fit exact recovery (2D graphene)" begin
        for tbm in (tbm_nn, tbm_nnn)
            Random.seed!(42)
            cs_r = randn(length(tbm))
            Em_r = spectrum(tbm(cs_r), ks)
            ptbm_fit = fit(tbm, Em_r, ks)
            @test ptbm_fit isa ParameterizedTightBindingModel
            @test rms(spectrum(ptbm_fit, ks) - Em_r) < 1e-5
        end
    end

    # -------------------------------------------------------------------------------------
    @testset "fit exact recovery (3D, SG 221)" begin
        # the canonical example from `fit`'s docstring, on random k-points (avoids a
        # Brillouin.jl dependency in the test environment)
        Random.seed!(123)
        cs_r = randn(length(tbm221))
        ks221 = [rand(3) for _ in 1:12]
        Em_r = spectrum(tbm221(cs_r), ks221)
        ptbm_fit = fit(tbm221, Em_r, ks221)
        @test rms(spectrum(ptbm_fit, ks221) - Em_r) < 1e-5
        @test ptbm_fit.cs ≈ cs_r          # amplitudes uniquely recovered for this model
    end

    # -------------------------------------------------------------------------------------
    @testset "rank-deficient moment seed (traceless terms)" begin
        # regression: terms 2 & 4 of `tbm221` are purely off-diagonal, so tr hᵢ(k) ≡ 0 and
        # the first-moment design matrix A[κ,i] = tr hᵢ(kₖ) has zero columns. The old seed
        # `A \ m₁` threw `SingularException` as soon as A was square (#k-points ≤ #terms);
        # the `pinv` seed is rank-robust and (correctly) seeds the traceless terms at 0.
        ktest = [rand(3) for _ in 1:4]
        for i in (2, 4)                       # confirm the precondition (traceless terms)
            @test all(k -> abs(tr(tbm221[i](k))) < 1e-12, ktest)
        end

        Random.seed!(221)
        cs_r  = randn(length(tbm221))
        ks_sq = [rand(3) for _ in 1:length(tbm221)]  # #k-points == #terms → square A
        Em_sq = spectrum(tbm221(cs_r), ks_sq)

        m = spectralmoments(tbm221, Em_sq, ks_sq)
        @test all(isfinite, m.c₀)             # threw `SingularException` before the pinv fix
        @test abs(m.c₀[2]) < 1e-10            # traceless term: no first-moment info → 0 seed
        @test abs(m.c₀[4]) < 1e-10

        # the full fit must also run (and recover) with this square, rank-deficient setup
        ptbm = fit(tbm221, Em_sq, ks_sq)
        @test rms(spectrum(ptbm, ks_sq) - Em_sq) < 1e-5
    end

    # -------------------------------------------------------------------------------------
    @testset "keyword arguments" begin
        Random.seed!(7)
        cs_r = randn(length(tbm_nnn))
        Em_r = spectrum(tbm_nnn(cs_r), ks)

        # `init`: seeding the first (deterministic) trial at the true answer converges
        # immediately, even with a single trial and no polish
        ptbm = fit(tbm_nn, spectrum(tbm_nn(cs_r[1:length(tbm_nn)]), ks), ks;
                   init = cs_r[1:length(tbm_nn)], max_multistarts = 1, polish = false)
        @test rms(spectrum(ptbm, ks) - spectrum(tbm_nn(cs_r[1:length(tbm_nn)]), ks)) < 1e-6

        # `lasso`: runs and returns a valid model (sparsity-encouraging penalty)
        ptbm_lasso = fit(tbm_nnn, Em_r, ks; lasso = 1e-3)
        @test ptbm_lasso isa ParameterizedTightBindingModel

        # `polish = false` and a first-order optimizer also fit (exercises non-default paths)
        ptbm_lbfgs = fit(tbm_nn, spectrum(tbm_nn(cs_r[1:length(tbm_nn)]), ks), ks;
                         optimizer = LBFGS(), polish = false)
        @test ptbm_lbfgs isa ParameterizedTightBindingModel
    end

    @testset "Broken example due to near-singular Hessian" begin
        # The following example errors due to an Optim.jl bug. I proposed a PR to fix the
        # bug at https://github.com/JuliaNLSolvers/Optim.jl/pull/1266, but in the meantime,
        # we just log it as a broken test

        # SG 221: (3d|A₁g) ⊕ (3d|B₂g); 6 bands, 3D (the `fit` docstring example)
        brs = calc_bandreps(221, Val(3))
        tbm = tb_hamiltonian(@composite brs[1] + brs[7])

        Random.seed!(0) # this value happens to trigger the bug; seed 1 does not e.g.
        ks = [rand(dim(brs)) for _ in 1:3]
        Em_ref = spectrum(tbm(randn(length(tbm))), ks)
        @test_broken fit(tbm, Em_ref, ks) isa ParameterizedTightBindingModel{3}
        # ↑ remove `@test_broken` once the Optim.jl#1266 PR is merged and released (& we 
        #   increase our compat)

        Random.seed!(1) # works on seed 1
        ks = [rand(dim(brs)) for _ in 1:3]
        Em_ref = spectrum(tbm(randn(length(tbm))), ks)
        @test fit(tbm, Em_ref, ks) isa ParameterizedTightBindingModel{3}
    end
end
