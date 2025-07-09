using SymmetricTightBinding, Test

include("pg_tb_hamiltonian.jl")    # plane groups
include("sg_tb_hamiltonian.jl")    # space groups
include("site_representations.jl") # site representations
include("symmetry-breaking.jl")    # symmetry breaking

#include("symmetry_analysis.jl")    # check that each tb model is symmetry compatible with
#                                    the constituents EBRs 

@testset "TightBindingModel" begin
    brs = calc_bandreps(16, Val(2))
    cbr = @composite brs[3]
    tbm = tb_hamiltonian(cbr, [[0,0],[1,0]])
    @testset "AbstractArray indexing into TightBindingModel" begin

        @test length(tbm) == 4
        tbm_subset1 = tbm[1:3]
        tbm_subset2 = tbm[[1,2,3]]
        @test length(tbm_subset1) == 3
        @test tbm_subset1 == tbm_subset2
        tbm_subset3 = tbm[[1, 4]]
        @test length(tbm_subset3) == 2
        @test tbm_subset3[1] == tbm[1] && tbm_subset3[2] == tbm[4]
    end

    @testset "Concatenation of TightBindingModels (`vcat`)" begin
        @test vcat(tbm[1:2], tbm[3:4]) == tbm
        @test vcat(tbm[1:2], tbm[3:3]) == tbm[1:3]
        @test vcat(tbm[1:2], tbm[4:4]) != tbm
    end
end
