using Test
using Crystalline
using SymmetricTightBinding
using SymmetricTightBinding: _GGRContribution
using StaticArrays: SVector
using LinearAlgebra: dot
using Random

# trapezoidal integration helper
_trapz(x, y) = sum(i -> (x[i+1] - x[i]) * (y[i] + y[i+1]) / 2, 1:length(x)-1)

@testset "Density of states" begin
    # ------------------------------------------------------------------------------------ #
    @testset "Per-cell GR contribution `_GGRContribution`" begin
        # Validate the closed-form per-cell contribution `C(δE)` against the distribution of
        # X = v·Δk for Δk uniform over the cube [-h,h]^D. Since Δk is uniform with density
        # 1/(2h)^D, the pdf of X is p(x) = C(x)/(2h)^D; we compare the analytic cumulative
        # F(t) = (2h)^{-D} ∫_{-wmax}^{t} C dδE against the empirical CDF of sampled X. This
        # independently pins down the (3D) coefficients and the meaning of |v| (the public API
        # only exercises `C` indirectly, so we test it directly here). We also check the mass
        # ∫C dδE = (2h)^D and the symmetry C(δE) = C(-δE).
        Random.seed!(1234)
        h = 0.0317
        for D in (1, 2, 3)
            for _ in 1:4
                v = ntuple(_ -> randn(), D)
                C = _GGRContribution(v, h)
                wmax = h * sum(abs, v)
                cellmass = (2h)^D

                # mass + symmetry (fine analytic grid)
                δs = range(-1.05wmax, 1.05wmax, 4001)
                Cs = [C(δ) for δ in δs]
                @test _trapz(δs, Cs) ≈ cellmass rtol=1e-3
                @test C(0.37wmax) ≈ C(-0.37wmax)

                # MC empirical CDF vs analytic CDF (Δk uniform in [-h,h]^D)
                N = 200_000
                vsv = SVector{D,Float64}(v)
                Xs = [dot(vsv, SVector{D,Float64}(ntuple(_ -> (2h)*(rand()-0.5), D)))
                      for _ in 1:N]
                gridδ = range(-wmax, wmax, 2001)
                gridC = [C(δ) for δ in gridδ]
                for t in range(-0.8wmax, 0.8wmax, 9)
                    mask = gridδ .≤ t
                    Fa = _trapz(gridδ[mask], gridC[mask]) / cellmass
                    Fe = count(≤(t), Xs) / N
                    @test Fa ≈ Fe atol=0.015
                end
            end
        end

        # explicit easy cases
        @test _GGRContribution((2.0,), 0.1)(0.0) ≈ 0.5          # 1D plateau 1/|v| inside |δ|≤|v|h
        @test _GGRContribution((2.0,), 0.1)(0.21) == 0.0        # 1D outside support (|v|h = 0.2)
        @test _GGRContribution((3.0, 0.0), 0.1)(0.0) ≈ 2*0.1/3  # 2D plateau 2h/a (b = 0)
        @test _GGRContribution((1e-8,), 0.1)(0.0) ≈ 1e8         # near-flat: integrable 1/|v| spike
        @test _GGRContribution((0.0,), 0.1)(0.0) == Inf         # perfectly flat: Inf at δ = 0
        @test _GGRContribution((0.0,), 0.1)(0.05) == 0.0        #                 0   elsewhere
    end

    # ------------------------------------------------------------------------------------ #
    @testset "1D analytic DOS (1-band cosine, SG 2)" begin
        # A 1-band 1D model has a pure-cosine dispersion E(k) = E₀ + W·cos(2πk), whose exact
        # per-unit-cell DOS is the van Hove form g(E) = 1/(π√(W² − (E−E₀)²)) for |E−E₀| < W
        # (note ∫g dE = 1, matching the single band). We test the full pipeline against it.
        brs = calc_bandreps(2, Val(1))
        cbr = @composite brs[1] # 1 band
        ptbm = tb_hamiltonian(cbr, [[0], [1]])([0.3, 0.5])
        Es = [only(spectrum_single_k(ptbm, [k])) for k in range(-0.5, 0.5, 4001)]
        E₀ = (maximum(Es) + minimum(Es)) / 2
        W = (maximum(Es) - minimum(Es)) / 2
        g_exact(E) = 1 / (π * sqrt(W^2 - (E - E₀)^2))

        Es = E₀ .+ range(-0.75W, 0.75W, 41) # interior points (away from van Hove band edges)
        g = densityofstates(ptbm, Es; Nk=4000)         # binned (M ≥ 2) path
        @test all(isapprox.(g, g_exact.(Es); rtol=3e-2))

        # single-energy (M == 1) path returns the pointwise GR density; also ≈ exact
        E₁ = E₀ + 0.3W
        @test only(densityofstates(ptbm, [E₁]; Nk=4000)) ≈ g_exact(E₁) rtol=3e-2
    end

    # ------------------------------------------------------------------------------------ #
    @testset "Sum rule ∫g dE ≈ N_bands (2D, 3D)" begin
        # the per-unit-cell normalization: integrating the DOS over all energies recovers
        # the number of bands. (The 1D case is covered analytically above.)
        function dos_integral(ptbm; Nk, margin=0.5, npts=1500)
            D = dim(CompositeBandRep(ptbm))
            ks = [SVector{D,Float64}(ntuple(_ -> rand()-0.5, D)) for _ in 1:400]
            Es = reduce(vcat, [spectrum_single_k(ptbm, k) for k in ks])
            Es = range(minimum(Es)-margin, maximum(Es)+margin, npts)
            return _trapz(Es, densityofstates(ptbm, Es; Nk))
        end

        @testset "2D (graphene, pg 17)" begin
            brs = calc_bandreps(17, Val(2))
            cbr = @composite brs[5] # (2b|A₁), 2 bands
            ptbm = tb_hamiltonian(cbr, [[0,0]])([0.0, 1.0]) # zero on-site, unit NN hopping
            @test dos_integral(ptbm; Nk=120) ≈ 2 rtol=2e-2
        end

        @testset "3D (SG 2)" begin
            brs = calc_bandreps(2, Val(3); timereversal=true)
            cbr = @composite brs[1] + brs[end] # 2 bands
            tbm = tb_hamiltonian(cbr, [[1,0,0],[0,1,0],[0,0,1]])
            Random.seed!(1)
            ptbm = tbm(randn(length(tbm)))
            @test dos_integral(ptbm; Nk=24) ≈ 2 rtol=4e-2
        end
    end

    # ------------------------------------------------------------------------------------ #
    @testset "Graphene DOS shape (Dirac dip)" begin
        # the graphene DOS vanishes linearly at the Dirac point (E = 0) and peaks at the van
        # Hove singularities (E=±2 here): a qualitative shape check against the global maximum
        brs = calc_bandreps(17, Val(2))
        cbr = @composite brs[5]
        ptbm = tb_hamiltonian(cbr, [[0,0]])([0.0, 1.0]) # bands span [-6, 6]

        Es = range(-6.5, 6.5, 601)
        g = densityofstates(ptbm, Es; Nk=150)
        iDirac = argmin(abs.(Es))           # E ≈ 0
        @test g[iDirac] < 0.05 * maximum(g) # deep dip at the Dirac point vs the van Hove peak
        @test all(≥(0), g)
    end

    # ------------------------------------------------------------------------------------ #
    @testset "`transform` keyword (chain rule on velocity)" begin
        # DOS in a transformed abscissa φ = transform(E) must satisfy g_φ(φ)·φ'(E) = g_E(E)
        # for a monotonic transform (conservation of states). Test with φ = E + a·E³.
        brs = calc_bandreps(17, Val(2))
        cbr = @composite brs[5]
        ptbm = tb_hamiltonian(cbr, [[0,0]])([0.0, 1.0])

        a = 0.05
        tf  = E -> E + a*E^3
        dtf = E -> 1 + 3a*E^2

        Es = range(-2.5, 2.5, 51)
        gE = densityofstates(ptbm, Es; Nk=160)
        gφ = densityofstates(ptbm, tf.(Es); Nk=160, transform=tf)
        # compare on interior points (avoid band-edge singularities where linearization is poor)
        for i in 10:42
            @test gφ[i]*dtf(Es[i]) ≈ gE[i] rtol=5e-2
        end

        # a `transform` returning a non-`Float64` `Real` is coerced internally (no MethodError)
        g32 = densityofstates(ptbm, Float32.(Es); Nk=40, transform = E -> Float32(E))
        @test all(isfinite, g32)
    end

    # ------------------------------------------------------------------------------------ #
    @testset "`bands` keyword (partial DOS)" begin
        # partial DOSs must partition the total: restricting to a set of bands drops exactly
        # those bands' contributions and nothing else.
        brs = calc_bandreps(2, Val(3); timereversal=true)
        cbr = @composite brs[1] + brs[end] # 2 bands
        tbm = tb_hamiltonian(cbr, [[1,0,0],[0,1,0],[0,0,1]])
        Random.seed!(1)
        ptbm = tbm(randn(length(tbm)))
        Nᵇ = tbm.N

        Es = range(-8, 8, 201)
        dE = step(Es)
        g_all = densityofstates(ptbm, Es; Nk=16)

        # each single-band partial DOS carries unit weight, and together they reproduce the
        # total DOS pointwise
        g_each = [densityofstates(ptbm, Es; Nk=16, bands=n:n) for n in 1:Nᵇ]
        for g in g_each
            @test sum(g)*dE ≈ 1 rtol=5e-2
            @test all(≥(0), g)
        end
        @test sum(g_each) ≈ g_all

        # passing every band is identical to omitting the keyword
        @test densityofstates(ptbm, Es; Nk=16, bands=1:Nᵇ) ≈ g_all

        # composes with `transform`: restricting commutes with transforming, since the chain
        # rule is applied per band
        tf = E -> E + 0.05E^3
        g_lo_t = densityofstates(ptbm, Es; Nk=16, bands=1:1, transform=tf)
        g_hi_t = densityofstates(ptbm, Es; Nk=16, bands=2:Nᵇ, transform=tf)
        @test g_lo_t .+ g_hi_t ≈ densityofstates(ptbm, Es; Nk=16, transform=tf)
    end

    # ------------------------------------------------------------------------------------ #
    @testset "errors & basic invariants" begin
        brs = calc_bandreps(17, Val(2))
        cbr = @composite brs[5]
        tbm = tb_hamiltonian(cbr, [[0,0]])
        ptbm = tbm([0.0, 1.0])

        # output length matches input length and is nonnegative
        Es = range(-3, 3, 123)
        g = densityofstates(ptbm, Es; Nk=40)
        @test length(g) == length(Es)
        @test all(≥(0), g)

        # unsorted or duplicated energies are rejected with an informative error
        @test_throws ErrorException densityofstates(ptbm, reverse(Es))
        @test_throws ErrorException densityofstates(ptbm, [-1.0, 0.0, 0.0, 1.0])

        # NONHERMITIAN models are rejected
        tbm_NH = tb_hamiltonian(cbr, [[0,0]], Val(NONHERMITIAN))
        ptbm_NH = tbm_NH(zeros(length(tbm_NH)))
        @test_throws ErrorException densityofstates(ptbm_NH, Es)
    end
end
