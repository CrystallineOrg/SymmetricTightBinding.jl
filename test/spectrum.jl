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
        es = spectrum_single_k(ptbm, k)
        @test length(es) == 2  # 2 bands
        @test issorted(es)     # eigenvalues sorted

        # compare against direct eigvals
        H = ptbm(k)
        @test H isa Hermitian
        @test es ≈ eigvals(H)
    end

    @testset "Multiple k-points" begin
        ks = [[0.0, 0.0], [0.5, 0.0], [1/3, 1/3], [0.0, 0.5]]
        Es = spectrum(ptbm, ks)
        @test size(Es) == (4, 2)    # (Nk, Nbands)
        @test all(issorted, eachrow(Es))  # each row sorted

        # verify consistency with single-k method
        for (i, k) in enumerate(ks)
            @test Es[i, :] ≈ spectrum_single_k(ptbm, k)
        end
    end

    @testset "Transform keyword" begin
        k = [0.0, 0.0]
        es = spectrum_single_k(ptbm, k)
        es_sq = spectrum_single_k(ptbm, k; transform=x->x^2)
        @test es_sq ≈ es.^2

        ks = [[0.0, 0.0], [0.5, 0.0]]
        Es = spectrum(ptbm, ks; transform=abs)
        @test all(Es .≥ 0)
    end

    @testset "Graphene Dirac point" begin
        # at K = (1/3, 1/3), graphene should have a degeneracy
        es_K = spectrum_single_k(ptbm, [1/3, 1/3])
        @test es_K[1] ≈ es_K[2] atol=1e-10

        # at Gamma, bands should be split
        es_Γ = spectrum_single_k(ptbm, [0.0, 0.0])
        @test abs(es_Γ[2] - es_Γ[1]) > 1e-1
    end

    @testset "D = 1 convenience method" begin
        brs_1D = calc_bandreps(2, Val(1))
        cbr_1D = @composite brs_1D[1]+brs_1D[3] # 2 bands
        tbm_1D = tb_hamiltonian(cbr_1D, [[0], [1]])
        ptbm_1D = tbm_1D([0.1, 0.2, -0.3, 0.4, -0.5])
        # for 1D models, we allow passing a vector of real numbers, each interpreted as
        # independent momentum points
        ks_numbers = [0.0, 0.25, 0.5, 0.75]     # ::Vector 
        ks_points  = [[k] for k in ks_numbers]
        @test spectrum(ptbm_1D, ks_numbers) ≈ spectrum(ptbm_1D, ks_points)

        ks_numbers_range = range(0, 1/2, 10)    # ::StepRange (an `AbstractVector`)
        ks_points_range  = [[k] for k in ks_numbers_range]
        @test spectrum(ptbm_1D, ks_numbers_range) ≈ spectrum(ptbm_1D, ks_points_range)

        # we should still throw an informative error if we try to do this with a D≠1 model
        @test_throws ErrorException spectrum(ptbm, ks_numbers)
    end
end
