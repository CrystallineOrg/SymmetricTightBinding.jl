using Test
using SymmetricTightBinding
using SymmetricTightBinding: sgrep_induced_by_siteir_excl_phase
using Crystalline

@testset "Site representations" begin
    @testset "SG #221" begin
        sgnum = 221
        brs = calc_bandreps(sgnum, Val(3))
        cbr = @composite brs[6]

        gens = generators(num(cbr), SpaceGroup{3})
        sgrep = sgrep_induced_by_siteir_excl_phase.(Ref(cbr), gens)

        @test length(sgrep) == length(gens) == 5
        @test gens == generators(sgnum, SpaceGroup{3})

        @testset "Handwritten representations" begin
            @test sgrep[1] == Complex[
                -1.0 0.0 0.0
                0.0 -1.0 0.0
                0.0 0.0 1.0
            ]
            @test sgrep[2] == Complex[
                -1.0 0.0 0.0
                0.0 1.0 0.0
                0.0 0.0 -1.0
            ]
            @test sgrep[3] == Complex[
                0.0 0.0 1.0
                1.0 0.0 0.0
                0.0 1.0 0.0
            ]
            @test sgrep[4] == Complex[
                0.0 1.0 0.0
                1.0 0.0 0.0
                0.0 0.0 -1.0
            ]
            @test sgrep[5] == Complex[
                -1.0 0.0 0.0
                0.0 -1.0 0.0
                0.0 0.0 -1.0
            ]
        end
    end # SG 221

    @testset "SG #224" begin
        sgnum = 224
        brs = calc_bandreps(sgnum, Val(3))
        cbr = @composite brs[13] + brs[19]

        gens = generators(num(cbr), SpaceGroup{3})
        sgrep = sgrep_induced_by_siteir_excl_phase.(Ref(cbr), gens)

        @test length(sgrep) == length(gens) == 5
        @test gens == generators(sgnum, SpaceGroup{3})

        @testset "Handwritten representations" begin
            @test sgrep[1] ≈ Complex[
                0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
                1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0
                0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0
                0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0
                0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0
                0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0
                0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0
            ]
            @test sgrep[2] ≈ Complex[
                0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0
                0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0
                1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
                0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0
                0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0
                0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0
                0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0
            ]
            @test sgrep[3] ≈ Complex[
                1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0
                0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0
                0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
                0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0
                0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0
                0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0
                0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0
            ]
            @test sgrep[4] ≈ Complex[
                 0.0 -1.0  0.0  0.0  0.0  0.0  0.0  0.0
                -1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                 0.0  0.0 -1.0  0.0  0.0  0.0  0.0  0.0
                 0.0  0.0  0.0 -1.0  0.0  0.0  0.0  0.0
                 0.0  0.0  0.0  0.0  0.0 -1.0  0.0  0.0
                 0.0  0.0  0.0  0.0 -1.0  0.0  0.0  0.0
                 0.0  0.0  0.0  0.0  0.0  0.0 -1.0  0.0
                 0.0  0.0  0.0  0.0  0.0  0.0  0.0 -1.0
            ]
            @test sgrep[5] ≈ Complex[
               -1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
                0.0 -1.0  0.0  0.0  0.0  0.0  0.0  0.0
                0.0  0.0 -1.0  0.0  0.0  0.0  0.0  0.0
                0.0  0.0  0.0 -1.0  0.0  0.0  0.0  0.0
                0.0  0.0  0.0  0.0 -1.0  0.0  0.0  0.0
                0.0  0.0  0.0  0.0  0.0 -1.0  0.0  0.0
                0.0  0.0  0.0  0.0  0.0  0.0 -1.0  0.0
                0.0  0.0  0.0  0.0  0.0  0.0  0.0 -1.0
            ]
        end
    end # SG 224

    @testset "Point Group #2 (-1)" begin
        sgnum = 2
        brs = calc_bandreps(sgnum, Val(1))
        cbr = @composite brs[2] + brs[3]

        gens = generators(num(cbr), SpaceGroup{1})
        sgrep = sgrep_induced_by_siteir_excl_phase.(Ref(cbr), gens)

        @test length(gens) == length(sgrep) == 1
        @test gens == generators(sgnum, SpaceGroup{1})

        @testset "Handwritten representation" begin
            @test sgrep[1] == Complex[-1.0-0.0im 0.0; 0.0 1.0]
        end
    end

    @testset "Graphene" begin
        sgnum = 17
        brs = calc_bandreps(sgnum, Val(2))
        cbr = @composite brs[5]

        gens = generators(num(cbr), SpaceGroup{2})
        sgrep = sgrep_induced_by_siteir_excl_phase.(Ref(cbr), gens)

        @test length(gens) == length(sgrep) == 3
        @test gens == generators(sgnum, SpaceGroup{2})

        @testset "Handwritten representation" begin
            @test sgrep[1] ≈ Complex[1.0 0.0; 0.0 1.0]
            @test sgrep[2] ≈ Complex[0.0 1.0; 1.0 0.0]
            @test sgrep[3] ≈ Complex[1.0 0.0; 0.0 1.0]
        end
    end

    @testset "Plane Group #10 (4)" begin
        sgnum = 10
        brs = calc_bandreps(sgnum, Val(2))
        cbr = @composite brs[1] + brs[end]

        gens = generators(num(cbr), SpaceGroup{2})
        sgrep = sgrep_induced_by_siteir_excl_phase.(Ref(cbr), gens)

        @test length(gens) == length(sgrep) == 2
        @test gens == generators(sgnum, SpaceGroup{2})

        @testset "Handwritten representation" begin
            @test(sgrep[1] ≈ Complex[
                1.0 0.0  0.0  0.0
                0.0 1.0  0.0  0.0
                0.0 0.0 -1.0  0.0
                0.0 0.0  0.0 -1.0
            ], atol=1e-13)
            @test_broken(sgrep[2] ≈ Complex[
                1.0 0.0 0.0    0.0
                0.0 1.0 0.0    0.0
                0.0 0.0 1.0im  0.0
                0.0 0.0 0.0   -1.0im
            ], atol=1e-13)
        end
    end
end # Site representations
