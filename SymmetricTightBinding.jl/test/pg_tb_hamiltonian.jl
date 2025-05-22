using Test
using SymmetricTightBinding
using Crystalline
using Crystalline: constant
using LinearAlgebra

@testset "TB examples in plane groups" begin

    # TODO: add more cases apart from Graphene

    @testset "Graphene" begin
        D = 2
        sgnum = 17
        brs = calc_bandreps(sgnum, Val(D))
        cbr = @composite brs[5]

        Rs = [[0, 0]] # TODO: maybe consider adding more hopping terms

        hop_orbits = obtain_symmetry_related_hoppings(Rs, brs[5], brs[5])

        tb_model = tb_hamiltonian(cbr, Rs)

        @testset "Hopping terms" begin
            @test length(hop_orbits) == 2 # 2 hopping terms, 1. onsite and 2. 1st nn

            @testset "Onsite" begin
                onsite = hop_orbits[1]
                @test isapprox(onsite.representative, RVec([0, 0]))
                @test length(onsite.orbit) == 1 # there should be only one distance orbit [0,0]
                @test length(onsite.hoppings[1]) == orbit(group(br)) # number of orbitals == number of onsite hoppings
            end

            @testset "1st NN" begin
                first_nn = hop_orbits[2]
                @test length(first_nn.orbit) == 2 * 3 # there are 3 1st nn and we consider separately
                #                                       going forward and backward hopping terms
                lattice_vecs = directbasis(sgnum, D)
                first_nn_dist =
                    norm(cartesianize(constant(first_nn.orbit[1]), lattice_vecs))
                @test first_nn_dist ≈ 1 / sqrt(3) # distance between 1st nn in Graphene

                for δᵢ in constant.(first_nn.orbit)
                    @test norm(δᵢ) == first_nn_dist # all 1st nn distances should be equal
                end
            end
        end # Hopping terms

        @testset "Hamiltonian" begin
            @test size(tb_model) == size(hop_orbits) # symmetry independent tb model == number of hopping orbits

            tb_onsite = tb_model[1]
            # dimensions of the hamiltonian should be equal to the number of orbitals
            @test size(tb_onsite) == (occupation(brs[5]), occupation(brs[5]))

            # TODO: finish this when the proper `tb_hamiltonian` structure is implemented
        end # Hamiltonian
    end # Graphene
end # plane groups
