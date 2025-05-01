using TETB, Test
using Crystalline

@testset "Site representations" begin
    @testset "SG #221" begin
        sgnum, D = 221, 3
        brs = calc_bandreps(sgnum, Val(D))
        coefs = zeros(Int, length(brs))
        coefs[6] = 1
        cbr = CompositeBandRep(coefs, brs)

        gens = generators(num(cbr), SpaceGroup{D})
        sgrep = sgrep_induced_by_siteir_excl_phase.(Ref(cbr), gens)

        @test length(sgrep) == length(gens) == 5
        @test gens == generators(sgnum, SpaceGroup{D})

        @testset "Handwritten representations" begin
            @test sgrep[1] == Complex[
                -1.0+0.0im 0.0+0.0im 0.0+0.0im
                0.0+0.0im -1.0+0.0im 0.0+0.0im
                0.0+0.0im 0.0+0.0im 1.0+0.0im
            ]
            @test sgrep[2] == Complex[
                -1.0+0.0im 0.0+0.0im 0.0+0.0im
                0.0+0.0im 1.0+0.0im 0.0+0.0im
                0.0+0.0im 0.0+0.0im -1.0+0.0im
            ]
            @test sgrep[3] == Complex[
                0.0+0.0im 0.0+0.0im 1.0+0.0im
                1.0+0.0im 0.0+0.0im 0.0+0.0im
                0.0+0.0im 1.0+0.0im 0.0+0.0im
            ]
            @test sgrep[4] == Complex[
                0.0+0.0im 1.0+0.0im 0.0+0.0im
                1.0+0.0im 0.0+0.0im 0.0+0.0im
                0.0+0.0im 0.0+0.0im -1.0+0.0im
            ]
            @test sgrep[5] == Complex[
                -1.0+0.0im 0.0+0.0im 0.0+0.0im
                0.0+0.0im -1.0+0.0im 0.0+0.0im
                0.0+0.0im 0.0+0.0im -1.0+0.0im
            ]
        end
    end # SG 221

    @testset "SG #224" begin
        sgnum, D = 224, 3
        brs = calc_bandreps(sgnum, Val(D))
        coefs = zeros(Int, length(brs))
        coefs[[13, 19]] .= 1
        cbr = CompositeBandRep(coefs, brs)

        gens = generators(num(cbr), SpaceGroup{D})
        sgrep = sgrep_induced_by_siteir_excl_phase.(Ref(cbr), gens)

        @test length(sgrep) == length(gens) == 5
        @test gens == generators(sgnum, SpaceGroup{D})

        @testset "Handwritten representations" begin
            @test sgrep[1] == Complex[
                0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im
                1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im
                0+0im 0+0im 0+0im 0+0im 1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im
            ]
            @test sgrep[2] == Complex[
                0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im 0+0im 0+0im 0+0im 0+0im
                1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im
                0+0im 0+0im 0+0im 0+0im 1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im
            ]
            @test sgrep[3] == Complex[
                1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0+0im 0+0im 0+0im 0+0im 1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im
            ]
            @test sgrep[4] == Complex[
                0.0+0.0im -1.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im
                -1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0.0+0.0im 0.0+0.0im -1.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0.0+0.0im 0.0+0.0im 0.0+0.0im -1.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im -1.0+0.0im 0.0+0.0im 0.0+0.0im
                0+0im 0+0im 0+0im 0+0im -1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im -1.0+0.0im 0.0+0.0im
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 0.0+0.0im -1.0+0.0im
            ]
            @test sgrep[5] == Complex[
                -1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0.0+0.0im -1.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0.0+0.0im 0.0+0.0im -1.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0.0+0.0im 0.0+0.0im 0.0+0.0im -1.0+0.0im 0+0im 0+0im 0+0im 0+0im
                0+0im 0+0im 0+0im 0+0im -1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im -1.0+0.0im 0.0+0.0im 0.0+0.0im
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im -1.0+0.0im 0.0+0.0im
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 0.0+0.0im -1.0+0.0im
            ]
        end
    end # SG 224

    @testset "Point Group #2 (-1)" begin
        sgnum, D = 2, 1
        brs = calc_bandreps(sgnum, Val(D))
        coefs = zeros(Int, length(brs))
        coefs[[2, 3]] .= 1
        cbr = CompositeBandRep(coefs, brs)

        gens = generators(num(cbr), SpaceGroup{D})
        sgrep = sgrep_induced_by_siteir_excl_phase.(Ref(cbr), gens)

        @test length(gens) == length(sgrep) == 1
        @test gens == generators(sgnum, SpaceGroup{D})

        @testset "Handwritten representation" begin
            @test sgrep[1] == Complex[-1.0-0.0im 0+0im; 0+0im 1.0+0.0im]
        end
    end

    @testset "Graphene" begin
        sgnum, D = 17, 2
        brs = calc_bandreps(sgnum, Val(D))
        coefs = zeros(Int, length(brs))
        coefs[5] = 1
        cbr = CompositeBandRep(coefs, brs)

        gens = generators(num(cbr), SpaceGroup{D})
        sgrep = sgrep_induced_by_siteir_excl_phase.(Ref(cbr), gens)

        @test length(gens) == length(sgrep) == 3
        @test gens == generators(sgnum, SpaceGroup{D})

        @testset "Handwritten representation" begin
            @test sgrep[1] == Complex[1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im]
            @test sgrep[2] == Complex[0.0+0.0im 1.0+0.0im; 1.0+0.0im 0.0+0.0im]
            @test sgrep[3] == Complex[1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im]
        end
    end

    @testset "Plane Group #10 (4)" begin
        sgnum = 10
        brs = calc_bandreps(sgnum, Val(2))
        coefs = zeros(Int, length(brs))
        coefs[[1, end]] .= 1
        cbr = CompositeBandRep(coefs, brs)

        gens = generators(num(cbr), SpaceGroup{D})
        sgrep = sgrep_induced_by_siteir_excl_phase.(Ref(cbr), gens)

        @test length(gens) == length(sgreps) == 2
        @test gens == generators(sgnum, SpaceGroup{D})

        @testset "Handwritten representation" begin
            @test sgreps[1] == Complex[
                1.0+0.0im 0.0+0.0im 0+0im 0+0im
                0.0+0.0im 1.0+0.0im 0+0im 0+0im
                0+0im 0+0im -1.0+0.0im 0.0+0.0im
                0+0im 0+0im 0.0+0.0im -1.0+0.0im
            ]
            @test sgreps[2] == Complex[
                0.0+0.0im 1.0+0.0im 0+0im 0+0im
                1.0+0.0im 0.0+0.0im 0+0im 0+0im
                0+0im 0+0im 0.0+1.0im 0.0+0.0im
                0+0im 0+0im 0.0+0.0im 0.0-1.0im
            ]
        end
    end
end # Site representations
