using Test
using SymmetricTightBinding
using Crystalline
using LinearAlgebra

@testset "Gradients" begin
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
        Es_ref = spectrum(ptbm, k)
        for i in 1:length(cs)
            cs₊ = copy(cs); cs₊[i] += ε
            cs₋ = copy(cs); cs₋[i] -= ε
            Es₊ = spectrum(tbm(cs₊), k)
            Es₋ = spectrum(tbm(cs₋), k)
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
end
