using Test
using SymmetricTightBinding
using Crystalline
using Crystalline: constant
using LinearAlgebra

@testset "TB examples in plane groups" begin
    @testset "Graphene" begin
        sgnum = 17
        brs = calc_bandreps(sgnum, Val(2))
        br = brs[5]
        cbr = @composite brs[5]

        hop_orbits = obtain_symmetry_related_hoppings([[0, 0]], br, br)
        tb_model = tb_hamiltonian(cbr, [[0, 0]])

        @testset "Hopping terms" begin
            @test length(hop_orbits) == 2 # 2 hopping terms, 1. onsite and 2. 1st nn

            @testset "Onsite" begin
                onsite = hop_orbits[1]
                @test isapprox(onsite.representative, RVec([0, 0]))
                @test length(onsite.orbit) == 1 # there should be only one distance orbit [0,0]
                @test length(onsite.hoppings[1]) == length(orbit(group(br))) # number of orbitals == number of onsite hoppings
            end

            @testset "1st NN" begin
                first_nn = hop_orbits[2]
                # check: there are three 1st nn & we consider forward + backward hopping
                # terms separately
                @test length(first_nn.orbit) == 2 * 3
              
                Rs = directbasis(sgnum, Val(2))
                first_nn_dist = norm(cartesianize(first_nn.orbit[1](), Rs))
                @test first_nn_dist ≈ 1 / sqrt(3) # distance between 1st nn in Graphene

                for δᵢ in first_nn.orbit
                    δᵢ_dist = norm(cartesianize(δᵢ(), Rs))
                    @test norm(δᵢ_dist) ≈ first_nn_dist # all 1st nn distances should be equal
                end
            end
        end # @testset "Hopping terms"

        @testset "Hamiltonian" begin
            @test size(tb_model) == size(hop_orbits) # symmetry independent tb model == number of hopping orbits

            tb_onsite = tb_model[1]
            # dimensions of the hamiltonian should be equal to the number of orbitals
            @test size(tb_onsite) == (occupation(br), occupation(br))

            # evaluate H(k) at a few k-points and check hermiticity
            ptbm = tb_model([0.0, 1.0])
            for k in [[0.0, 0.0], [0.1, 0.2], [1/3, 1/3]]
                H = ptbm(k)
                @test H ≈ H'
            end
        end # @testset "Hamiltonian"

    end # @testset "Graphene"

    @testset "Square lattice (PG 10, p4mm)" begin
        sgnum = 11 # p4mm
        brs = calc_bandreps(sgnum, Val(2))
        cbr = @composite brs[1]
        tbm = tb_hamiltonian(cbr, [[0, 0]])

        cs = [0.3*cospi(0.73*n) for n in 1:length(tbm)] # deterministic "random" coefficients
        ptbm = tbm(cs)
        # C₄ symmetry: H(k₁,k₂) and H(k₂,-k₁) should have same spectrum
        k = [0.1, 0.3]
        k_rot = [0.3, -0.1] # C₄ rotation in reciprocal space
        es = spectrum_single_k(ptbm, k)
        es_rot = spectrum_single_k(ptbm, k_rot)
        @test es ≈ es_rot  atol=1e-10
    end

    @testset "Triangular lattice (PG 13, p3)" begin
        sgnum = 13 # p3
        brs = calc_bandreps(sgnum, Val(2))
        cbr = @composite brs[1]
        tbm = tb_hamiltonian(cbr, [[0, 0]])
        @test length(tbm) ≥ 1

        if length(tbm) > 0
            cs = [0.3*cospi(0.73*n) for n in 1:length(tbm)]
            ptbm = tbm(cs)
            H = ptbm([-0.1, 0.3])
            @test H ≈ H'
        end
    end

    @testset "Oblique lattice (PG 2, p2)" begin
        sgnum = 2 # p2
        brs = calc_bandreps(sgnum, Val(2))
        cbr = @composite brs[1]
        tbm = tb_hamiltonian(cbr, [[0, 0], [1, 0]])
        @test length(tbm) > 0

        cs = [0.3*cospi(0.73*n) for n in 1:length(tbm)]
        ptbm = tbm(cs)
        # TRS: H(k) = H*(-k)
        k = [0.1, 0.2]
        H_k₊ = copy(ptbm(k))
        H_k₋ = copy(ptbm(-k))
        @test H_k₊ ≈ conj.(H_k₋)
    end

end # @testset "TB examples in plane groups"
