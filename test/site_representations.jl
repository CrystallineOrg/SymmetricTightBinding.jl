using Test
using SymmetricTightBinding
using SymmetricTightBinding: site_induced_sgrep_excl_phase, site_induced_timereversal_unitary,
                             OrbitalOrdering, obtain_basis_free_parameters
using Crystalline
using LinearAlgebra: I

@testset "Site representations" begin
    @testset "SG #221" begin
        sgnum = 221
        brs = bandreps(sgnum, Val(3))
        cbr = @composite brs[6]

        gens = generators(num(cbr), SpaceGroup{3})
        sgrep = site_induced_sgrep_excl_phase.(Ref(cbr), gens)

        @test length(sgrep) == length(gens) == 5
        @test gens == generators(sgnum, SpaceGroup{3})

        @testset "Handwritten representations" begin
            @test sgrep[1] == Complex[
               -1.0  0.0 0.0
                0.0 -1.0 0.0
                0.0  0.0 1.0
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
                0.0 1.0  0.0
                1.0 0.0  0.0
                0.0 0.0 -1.0
            ]
            @test sgrep[5] == Complex[
               -1.0  0.0  0.0
                0.0 -1.0  0.0
                0.0  0.0 -1.0
            ]
        end
    end # SG 221

    @testset "SG #224" begin
        sgnum = 224
        brs = bandreps(sgnum, Val(3))
        cbr = @composite brs[13] + brs[19]

        gens = generators(num(cbr), SpaceGroup{3})
        sgrep = site_induced_sgrep_excl_phase.(Ref(cbr), gens)

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
        brs = bandreps(sgnum, Val(1))
        cbr = @composite brs[2] + brs[3]

        gens = generators(num(cbr), SpaceGroup{1})
        sgrep = site_induced_sgrep_excl_phase.(Ref(cbr), gens)

        @test length(gens) == length(sgrep) == 1
        @test gens == generators(sgnum, SpaceGroup{1})

        @testset "Handwritten representation" begin
            @test sgrep[1] == Complex[-1.0-0.0im 0.0; 0.0 1.0]
        end
    end

    @testset "Graphene" begin
        sgnum = 17
        brs = bandreps(sgnum, Val(2))
        cbr = @composite brs[5]

        gens = generators(num(cbr), SpaceGroup{2})
        sgrep = site_induced_sgrep_excl_phase.(Ref(cbr), gens)

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
        brs = bandreps(sgnum, Val(2))
        cbr = @composite brs[1] + brs[end]

        gens = generators(num(cbr), SpaceGroup{2})
        sgrep = site_induced_sgrep_excl_phase.(Ref(cbr), gens)

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

    @testset "Orbital ordering (#141)" begin
        # (3g|Egˢ) in P6/mmm: 3 sites with 2 partner functions (a Kramers pair) each
        brs = bandreps(191; spinful = Val(true))
        br = brs[1]
        V, Q = length(orbit(group(br))), Crystalline.irdim(br.siteir)
        @test (V, Q) == (3, 2)

        # default ordering is partner-function-major, i.e., sites run fastest
        ordering = OrbitalOrdering(br)
        @test [o.site_idx for o in ordering] == repeat(1:V, Q)
        @test [o.partner_idx for o in ordering] == repeat(1:Q; inner = V)
        # positions follow the ordering (i.e., agree with its `wp`s, as used in `Mm`)
        @test orbital_positions(br) ≈ [constant(parent(o.wp)) for o in ordering]
        Γ = site_induced_timereversal_unitary(br)
        @test Γ == kron(timereversal_unitary(br.siteir), I(V)) # = iσʸ ⊗ 𝟙, pseudospin-blocked

        # every ordering-dependent quantity must permute consistently with the ordering: check
        # this for a site-major ordering
        p = [something(findfirst(o -> o.site_idx == i && o.partner_idx == j, ordering)) 
                                                                for i in 1:V for j in 1:Q]
        ordering′ = OrbitalOrdering(ordering.ordering[p])
        @test [o.site_idx for o in ordering′] == repeat(1:V; inner = Q)
        @test site_induced_timereversal_unitary(br, ordering′) == Γ[p, p]
        for g in SymmetricTightBinding.primitivized_generators(br)
            ρ = site_induced_sgrep_excl_phase(br, g)
            @test site_induced_sgrep_excl_phase(br, g, ordering′) == ρ[p, p]
        end

        # the coefficient basis is indexed by partner functions, not orbitals, so it does not
        # depend on the ordering; the M-matrix is simply permuted
        h_orbits = SymmetricTightBinding.obtain_symmetry_related_hoppings(
            [[0, 0, 0], [1, 0, 0]], br, br)
        for h_orbit in h_orbits
            Mm, ts = obtain_basis_free_parameters(h_orbit, br, br)
            Mm′, ts′ = obtain_basis_free_parameters(h_orbit, br, br, ordering′, ordering′)
            @test Mm′ == Mm[:, :, p, p]
            @test ts′ ≈ ts
        end
    end
end # Site representations
