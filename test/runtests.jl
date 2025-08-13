using SymmetricTightBinding, Test

include("pg_tb_hamiltonian.jl")    # plane groups
include("sg_tb_hamiltonian.jl")    # space groups
include("site_representations.jl") # site representations
include("symmetry-breaking.jl")    # symmetry breaking

#include("symmetry_analysis.jl")    # check that each tb model is symmetry compatible with
#                                    the constituents EBRs 

@testset "AbstractArray interface" begin
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

@testset "Issue #73: multi-EBR without TR" begin
    sgnum = 22
    brs = calc_bandreps(sgnum, Val(3); timereversal=false)
    cbr = @composite brs[9]+brs[9+4] # (4b|A) + (4a|A) (2 bands)
    # simply test that it doesn't error
    @test tb_hamiltonian(cbr, [[1,1,1]]) isa TightBindingModel
end

@testset "Issue #85" begin
    sgnum = 4
    brs = calc_bandreps(sgnum, Val(2))

    # Sub-issue 1: Wyckoff positions in orbit outside primitive unit cell [0,1)ᴰ
    br_small_αβγ = SymmetricTightBinding.pin_free(brs[1], [.1, .2])
    orbit_small_αβγ = orbit(group(br_small_αβγ))
    @test all(q -> all(qᵢ -> 0 ≤ qᵢ < 1, q()), orbit_small_αβγ)

    br_big_αβγ = SymmetricTightBinding.pin_free(brs[1], [.9, .8])
    orbit_big_αβγ = orbit(group(br_big_αβγ))
    @test all(q -> all(qᵢ -> 0 ≤ qᵢ < 1, q()), orbit(group(br_big_αβγ)))

    # Sub-issue 2: missing merging of time-reversal related "partial" orbits
    for br in [br_small_αβγ, br_big_αβγ]
        h_orbits = obtain_symmetry_related_hoppings([[0,0]], br, br)
        @test length(h_orbits) == 2
        @test length(h_orbits[1].orbit) == 1 # on-site
        @test length(h_orbits[2].orbit) == 4 # hopping; combination of two orbits, merged by TR
    end
end