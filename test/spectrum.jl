using Test
using SymmetricTightBinding
using Crystalline
using LinearAlgebra

@testset "Spectrum" begin
    # set up graphene model (plane group 17, (2b|A₁) EBR)
    brs = calc_bandreps(17, Val(2))
    cbr = @composite brs[5]
    tbm = tb_hamiltonian(cbr, [[0, 0]])
    cs = [0.0, 1.0]
    ptbm = tbm(cs)

    @testset "Single k-point" begin
        k = [0.2, -0.1]
        es = spectrum(ptbm, k)
        @test length(es) == 2  # 2 bands
        @test issorted(es)     # eigenvalues sorted

        # compare against direct eigvals
        H = Hermitian(ptbm(k))
        @test es ≈ eigvals(H)
    end

    @testset "Multiple k-points" begin
        ks = [[0.0, 0.0], [0.5, 0.0], [1/3, 1/3], [0.0, 0.5]]
        Es = spectrum(ptbm, ks)
        @test size(Es) == (4, 2)    # (Nk, Nbands)
        @test all(issorted, eachrow(Es))  # each row sorted

        # verify consistency with single-k method
        for (i, k) in enumerate(ks)
            @test Es[i, :] ≈ spectrum(ptbm, k)
        end
    end

    @testset "Transform keyword" begin
        k = [0.0, 0.0]
        es = spectrum(ptbm, k)
        es_sq = spectrum(ptbm, k; transform=x->x^2)
        @test es_sq ≈ es.^2

        ks = [[0.0, 0.0], [0.5, 0.0]]
        Es = spectrum(ptbm, ks; transform=abs)
        @test all(Es .≥ 0)
    end

    @testset "Graphene Dirac point" begin
        # at K = (1/3, 1/3), graphene should have a degeneracy
        es_K = spectrum(ptbm, [1/3, 1/3])
        @test es_K[1] ≈ es_K[2] atol=1e-10

        # at Gamma, bands should be split
        es_Γ = spectrum(ptbm, [0.0, 0.0])
        @test abs(es_Γ[2] - es_Γ[1]) > 1e-1
    end
end
