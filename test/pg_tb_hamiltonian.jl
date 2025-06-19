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

            # TODO: finish this when the proper `tb_hamiltonian` structure is implemented
        end # @testset "Hamiltonian"

    end # @testset "Graphene"

    # TODO: add more cases apart from Graphene

end # @testset "TB examples in plane groups"
