using Test
using SymmetricTightBinding
using SymmetricTightBinding: _degenerate_band_groups
using Crystalline
using LinearAlgebra

@testset "Gradients" begin
    @testset "_degenerate_band_groups" begin
        # `Es` is assumed sorted ascending (as `solve` returns). With well-separated energies
        # all bands are singletons; (near-)degenerate runs of any multiplicity must group into
        # a single range — in particular 3-fold and higher, not just pairs (regression test).
        groups(Es; rtol=1e-12, atol=1e-9) = _degenerate_band_groups(Es, rtol, atol)

        @test groups([1.0, 2.0, 3.0]) == [1:1, 2:2, 3:3]      # all non-degenerate
        @test groups([1.0, 1.0, 2.0]) == [1:2, 3:3]           # leading pair
        @test groups([1.0, 2.0, 2.0]) == [1:1, 2:3]           # trailing pair (last-band path)
        @test groups([1.0, 1.0, 1.0]) == [1:3]                # 3-fold (previously split as 1:2,3:3)
        @test groups([1.0, 1.0, 1.0, 1.0]) == [1:4]           # 4-fold
        @test groups([0.0, 1.0, 1.0, 1.0, 2.0]) == [1:1, 2:4, 5:5] # 3-fold in the middle
        @test groups([0.0, 1.0, 1.0, 1.0]) == [1:1, 2:4]      # 3-fold ending the spectrum
        @test groups([1.0]) == [1:1]                          # single band

        # near-degeneracy within `atol` is grouped; gaps above it are not
        @test groups([1.0, 1.0 + 1e-13, 1.0 + 2e-13]) == [1:3]
        @test groups([1.0, 1.0 + 1e-6, 1.0 + 2e-6]; atol=1e-9) == [1:1, 2:2, 3:3]
    end

    # set up graphene model
    brs = calc_bandreps(17, Val(2))
    cbr = @composite brs[5]
    tbm = tb_hamiltonian(cbr, [[0, 0]])
    cs = [0.0, 1.0]
    ptbm = tbm(cs)

    @testset "gradient_wrt_hopping" begin
        tbmg = gradient_wrt_hopping(tbm)

        k = [0.1, 0.2]
        ∇H = tbmg(k)
        @test length(∇H) == length(tbm) # one matrix per term/coefficient

        # each gradient component should be an N×N matrix
        N = tbm.N
        for dH in ∇H
            @test size(dH) == (N, N)
        end

        # verify: H(k) = ∑ cᵢ ∂H(k)/∂cᵢ (cf. H(k) = ∑ cᵢ Hᵢ(k))
        H_reconstructed = sum(c * dH for (c, dH) in zip(cs, ∇H))
        H_direct = ptbm(k)
        @test H_reconstructed ≈ H_direct

        # single-index access
        dH1 = tbmg(k, 1)
        @test dH1 ≈ ∇H[1]

        # gradient from ptbm should equal gradient from tbm
        tbmg2 = gradient_wrt_hopping(ptbm)
        @test tbmg2(k) ≈ ∇H  atol=1e-14
    end

    @testset "gradient_wrt_momentum" begin
        ∇ptbm = gradient_wrt_momentum(ptbm)

        k = [0.1, 0.2]
        ∇Hs = ∇ptbm(k, (1, 2))
        @test length(∇Hs) == 2  # D=2 components

        N = tbm.N
        for dH in ∇Hs
            @test size(dH) == (N, N)
        end

        # verify via finite differences: ∂H/∂kᵢ ≈ [H(k+εeᵢ) - H(k-εeᵢ)] / 2ε
        ε = 1e-7
        for i in 1:2
            dk = zeros(2)
            dk[i] = ε
            H₊  = copy(ptbm(k .+ dk))  # copy: ptbm returns mutable scratch
            H₋ = copy(ptbm(k .- dk))
            dH_fd = (H₊ .- H₋) ./ (2ε)
            @test ∇Hs[i] ≈ dH_fd  atol=1e-4
        end
    end

    @testset "energy_gradient_wrt_hopping" begin
        k = [0.1, 0.2]
        ∇Es = energy_gradient_wrt_hopping(ptbm, k)
        @test length(∇Es) == tbm.N  # one gradient per band

        # each gradient vector has length = number of coefficients
        for ∇E in ∇Es
            @test length(∇E) == length(tbm)
        end

        # verify via finite differences: ∂Eₙ/∂cᵢ ≈ [Eₙ(c+εeᵢ) - Eₙ(c-εeᵢ)] / 2ε
        ε = 1e-8
        Es_ref = spectrum_single_k(ptbm, k)
        for i in 1:length(cs)
            cs₊ = copy(cs); cs₊[i] += ε
            cs₋ = copy(cs); cs₋[i] -= ε
            Es₊ = spectrum_single_k(tbm(cs₊), k)
            Es₋ = spectrum_single_k(tbm(cs₋), k)
            dEs_fd = (Es₊ .- Es₋) ./ (2ε)
            dEs_analytic = [∇Es[n][i] for n in 1:tbm.N]
            @test dEs_analytic ≈ dEs_fd  atol=1e-5
        end
    end

    @testset "Degenerate energy gradient" begin
        # at the Dirac point K=(1/3,1/3), graphene has a degeneracy
        k_K = [1/3, 1/3]
        ∇Es = energy_gradient_wrt_hopping(ptbm, k_K)
        @test length(∇Es) == tbm.N

        # test against known value of gradient [1, 0] here, for both bands, cf. degeneracy
        @test ∇Es ≈ [[1.0, 0.0], [1.0, 0.0]]  atol=1e-14
    end

    @testset "energy_gradient_wrt_momentum" begin
        k = [0.1, 0.2]
        vs = energy_gradient_wrt_momentum(ptbm, k)
        @test length(vs) == tbm.N          # one group velocity per band
        @test all(v -> length(v) == 2, vs) # each velocity is a D=2 vector

        # verify via finite differences: ∂Eₙ/∂kᵢ ≈ [Eₙ(k+εeᵢ) - Eₙ(k-εeᵢ)] / 2ε
        ε = 1e-7
        for i in 1:2
            dk = zeros(2); dk[i] = ε
            Es₊ = spectrum_single_k(ptbm, k .+ dk)
            Es₋ = spectrum_single_k(ptbm, k .- dk)
            dEs_fd = (Es₊ .- Es₋) ./ (2ε)
            dEs_analytic = [vs[n][i] for n in 1:tbm.N]
            @test dEs_analytic ≈ dEs_fd  atol=1e-5
        end
    end

    @testset "Degenerate momentum gradient (Dirac point)" begin
        # at the Dirac point K=(1/3,1/3) graphene is degenerate; the degenerate variant must
        # still return one finite velocity per band. Diagonalizing each momentum component in
        # the 2D Dirac subspace yields opposite-signed eigenvalues, so the two bands' velocities
        # come out as exact negatives of one another (the two halves of the linear crossing).
        k_K = [1/3, 1/3]
        vs = energy_gradient_wrt_momentum(ptbm, k_K)
        @test length(vs) == tbm.N
        @test all(v -> all(isfinite, v), vs)
        @test vs[1] ≈ -vs[2]  atol=1e-12
    end
end
