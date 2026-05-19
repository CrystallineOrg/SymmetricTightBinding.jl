using Test
using SymmetricTightBinding
using Crystalline

@testset "Symmetry breaking" begin
    @testset "2D example from docs" begin
        D = 2
        brs = calc_bandreps(11, Val(D); timereversal = true)
        cbr = @composite brs[1] # (2c|A₁)
        tbm = tb_hamiltonian(cbr, [[0,0], [1,0]])
        
        Δtbm_m   = subduced_complement(tbm, 10)                       # break mirror
        Δtbm_tr  = subduced_complement(tbm, 11; timereversal = false) # break TR
        Δtbm_mtr = subduced_complement(tbm, 10; timereversal = false) # break both
        @test length(Δtbm_m) == 1
        @test length(Δtbm_tr) == 1
        @test length(Δtbm_mtr) == 4
        
        # restrict to simpler terms (tbm[5] = complicated diagonally-directed hopping term)
        tbm′ = tbm[1:4] 
        Δtbm′_m   = subduced_complement(tbm′, 10)                       # break mirror
        Δtbm′_tr  = subduced_complement(tbm′, 11; timereversal = false) # break TR
        Δtbm′_mtr = subduced_complement(tbm′, 10; timereversal = false) # break both
        @test length(Δtbm′_m) == 0
        @test length(Δtbm′_tr) == 0
        @test length(Δtbm′_mtr) == 1

        Δtbm′_C4 = subduced_complement(tbm′, 6)               # break C₄
        @test length(Δtbm′_C4) == 3
        Δtbm_diagonal_term = subduced_complement(tbm[5:5], 6) # break C₄ (diagonal term only)
        @test length(Δtbm_diagonal_term) == 1
        
        # adding more orbits, we find a simple mirror-breaking term, but not a TR-breaking
        tbm_big_full = tb_hamiltonian(cbr, [[0,0], [1,0], [1,1]])
        tbm_big = tbm_big_full[[1:4...,6]] # drop complicated diagonally-directed term again
        Δtbm_big_m   = subduced_complement(tbm_big, 10)                       # break mirror
        Δtbm_big_tr  = subduced_complement(tbm_big, 11; timereversal = false) # break TR
        Δtbm_big_mtr = subduced_complement(tbm_big, 10; timereversal = false) # break both
        @test length(Δtbm_big_m) == 1
        @test length(Δtbm_big_tr) == 0
        @test length(Δtbm_big_mtr) == 2
        @test issubset(Δtbm_big_m, Δtbm_big_mtr) # must subset eachother

        # breaking mirror and TR symmetry together should give the same basis as starting
        # directly with plane group p4 (#10) and breaking TR from the get-go
        brs10 = calc_bandreps(10, Val(D); timereversal=false)
        cbr10 = @composite brs10[1] # (2c|A) (unlike #11: two distinct M irreps, M₃ & M₄)

        tbm10 = tb_hamiltonian(cbr10, [[0,0], [1,0], [0,1]]) # (⋆)
        @test length(tbm10) == length(tbm) + length(Δtbm_mtr)
        # (⋆): must add [0,1] also, cf. diagonally-directed hopping term involving both 
        # [1,0] and [0,1] in C₄ setting

        tbm10_big = tb_hamiltonian(cbr10, [[0,0], [1,0], [1,1]])
        Δtbm_big_full_mtr = subduced_complement(tbm_big_full, 10; timereversal = false)
        @test length(tbm10_big) == length(tbm_big_full) + length(Δtbm_big_full_mtr)
    end

    @testset "3D example" begin
        # this is not a well-thought out example, but just added to test that things work
        # without erroring
        brs = calc_bandreps(16, Val(3); timereversal = true)
        cbr = @composite brs[1] + brs[end] #  (1h|A) + (1a|B₂) (2 bands)
        tbm = tb_hamiltonian(cbr, [[0,0,0],])

        @test length(subduced_complement(tbm, 3)) == 1
    end

    @testset "broken 3D example" begin
        brs = calc_bandreps(121, Val(3); timereversal = true)
        cbr = @composite brs[1] + brs[end-1] # (4d|A) + (2a|A₂) (3 bands)
        tbm = tb_hamiltonian(cbr, [[0,0,0],])

        @test_broken subduced_complement(tbm, 82)
        # ERROR: failed to find any nonzero block (br=(4d|A), siteg=[1, {2₁₁₀|1,0,1},
        #        {-4₁₁₀⁺|1,0,0}, {-4₁₁₀⁻|1,1,1}] (4d: [3/4, 1/4, 1/2]), op=2₀₀₁)
    end
end