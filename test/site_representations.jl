using TETB, Test
using Crystalline

@testset "Site representations" begin

    @testset "SG #221" begin
        sg_num, D = 221, 3
        brs = calc_bandreps(sg_num, Val(D))
        coefs = zeros(Int, length(brs))
        coefs[6] = 1
        cbr = CompositeBandRep(coefs, brs)

        gens, sg_rep = sgrep_induced_by_siteir_generators(cbr)

        @test length(sg_rep) == length(gens) == 5
        @test gens == generators(sg_num, SpaceGroup{D})

        @testset "Handwritten representations" begin
            @test sg_rep[1] == Complex[-1.0+0.0im 0.0+0.0im 0.0+0.0im;
                0.0+0.0im -1.0+0.0im 0.0+0.0im;
                0.0+0.0im 0.0+0.0im 1.0+0.0im]
            @test sg_rep[2] == Complex[-1.0+0.0im 0.0+0.0im 0.0+0.0im;
                0.0+0.0im 1.0+0.0im 0.0+0.0im;
                0.0+0.0im 0.0+0.0im -1.0+0.0im]
            @test sg_rep[3] == Complex[0.0+0.0im 0.0+0.0im 1.0+0.0im;
                1.0+0.0im 0.0+0.0im 0.0+0.0im;
                0.0+0.0im 1.0+0.0im 0.0+0.0im]
            @test sg_rep[4] == Complex[0.0+0.0im 1.0+0.0im 0.0+0.0im;
                1.0+0.0im 0.0+0.0im 0.0+0.0im;
                0.0+0.0im 0.0+0.0im -1.0+0.0im]
            @test sg_rep[5] == Complex[-1.0+0.0im 0.0+0.0im 0.0+0.0im;
                0.0+0.0im -1.0+0.0im 0.0+0.0im;
                0.0+0.0im 0.0+0.0im -1.0+0.0im]
        end
    end # SG 221

    @testset "SG #224" begin
        sg_num, D = 224, 3
        brs = calc_bandreps(sg_num, Val(D))
        coefs = zeros(Int, length(brs))
        coefs[[13, 19]] .= 1
        cbr = CompositeBandRep(coefs, brs)

        gens, sg_rep = sgrep_induced_by_siteir_generators(cbr)

        @test length(sg_rep) == length(gens) == 5
        @test gens == generators(sg_num, SpaceGroup{D})

        @testset "Handwritten representations" begin
            @test sg_rep[1] == Complex[0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im;
                0+0im 0+0im 0+0im 0+0im 1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im;
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im;
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im]
            @test sg_rep[2] == Complex[0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im;
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im;
                0+0im 0+0im 0+0im 0+0im 1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im;
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im]
            @test sg_rep[3] == Complex[1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0+0im 0+0im 0+0im 0+0im 1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im;
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im;
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im;
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im]
            @test sg_rep[4] == Complex[0.0+0.0im -1.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                -1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0.0+0.0im 0.0+0.0im -1.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0.0+0.0im 0.0+0.0im 0.0+0.0im -1.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im -1.0+0.0im 0.0+0.0im 0.0+0.0im;
                0+0im 0+0im 0+0im 0+0im -1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im;
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im -1.0+0.0im 0.0+0.0im;
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 0.0+0.0im -1.0+0.0im]
            @test sg_rep[5] == Complex[-1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0.0+0.0im -1.0+0.0im 0.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0.0+0.0im 0.0+0.0im -1.0+0.0im 0.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0.0+0.0im 0.0+0.0im 0.0+0.0im -1.0+0.0im 0+0im 0+0im 0+0im 0+0im;
                0+0im 0+0im 0+0im 0+0im -1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im;
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im -1.0+0.0im 0.0+0.0im 0.0+0.0im;
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im -1.0+0.0im 0.0+0.0im;
                0+0im 0+0im 0+0im 0+0im 0.0+0.0im 0.0+0.0im 0.0+0.0im -1.0+0.0im]
        end
    end # SG 224

    @testset "Point Group #2 (-1)" begin
        pg_num, D = 2, 1
        brs = calc_bandreps(pg_num, Val(D))
        coefs = zeros(Int, length(brs))
        coefs[[2, 3]] .= 1
        cbr = CompositeBandRep(coefs, brs)

        gens, sg_rep = sgrep_induced_by_siteir_generators(cbr)

        @test length(gens) == length(sg_rep) == 1
        @test gens == generators(pg_num, SpaceGroup{D})

        @testset "Handwritten representation" begin
            @test sg_rep[1] == Complex[-1.0-0.0im 0+0im; 0+0im 1.0+0.0im]
        end
    end

    @testset "Graphene" begin
        pg_num, D = 17, 2
        brs = calc_bandreps(pg_num, Val(D))
        coefs = zeros(Int, length(brs))
        coefs[5] = 1
        cbr = CompositeBandRep(coefs, brs)

        gens, sg_rep = sgrep_induced_by_siteir_generators(cbr)

        @test length(gens) == length(sg_rep) == 3
        @test gens == generators(pg_num, SpaceGroup{D})

        @testset "Handwritten representation" begin
            @test sg_rep[1] == Complex[1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im]
            @test sg_rep[2] == Complex[0.0+0.0im 1.0+0.0im; 1.0+0.0im 0.0+0.0im]
            @test sg_rep[3] == Complex[1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im]
        end
    end

    @testset "Plane Group #10 (4)" begin
        pg_num = 10
        brs = calc_bandreps(pg_num, Val(2))
        coefs = zeros(Int, length(brs))
        coefs[[1, end]] .= 1
        cbr = CompositeBandRep(coefs, brs)

        gens, sg_reps = sgrep_induced_by_siteir_generators(cbr)

        @test length(gens) == length(sg_reps) == 2
        @test gens == generators(pg_num, SpaceGroup{D})

        @testset "Handwritten representation" begin
            @test sg_reps[1] == Complex[1.0+0.0im 0.0+0.0im 0+0im 0+0im;
                0.0+0.0im 1.0+0.0im 0+0im 0+0im;
                0+0im 0+0im -1.0+0.0im 0.0+0.0im;
                0+0im 0+0im 0.0+0.0im -1.0+0.0im]
            @test sg_reps[2] == Complex[0.0+0.0im 1.0+0.0im 0+0im 0+0im;
                1.0+0.0im 0.0+0.0im 0+0im 0+0im;
                0+0im 0+0im 0.0+1.0im 0.0+0.0im;
                0+0im 0+0im 0.0+0.0im 0.0-1.0im]
        end
    end
end # Site representations