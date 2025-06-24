using Test
using SymmetricTightBinding, Crystalline

@testset "Symmetry analysis" begin
    for D in 1:3
        for sgnum in MAX_SGNUM[D]
            @testset "Space group $sgnum in dimension $D" begin
                brs = calc_bandreps(sgnum, Val(D))
                for br in brs
                    cbr = @composite br
                    ptbm = tb_hamiltonian(cbr)(randn(length(cbr)))

                    v = SymmetryVector(cbr)
                    v_model = collect_compatible(ptbm)

                    @test v == v_model[1] # check that the symmetry vector is compatible with the model
                end # for br in brs
            end # @testset "Space group $sgnum in dimension $D"
        end # for sgnum in MAX_SGNUM[D]
    end # for D in 1:3
end # @testset "Symmetry analysis"