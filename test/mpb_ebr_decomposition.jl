using TETB, Test
using Crystalline # for calc_bandreps

@testset "EBR decomposition" begin

    @testset "SG #221" begin

        # construct the structure under study
        R1 = 0.2 # cylinder radius
        mat = mp.Medium(epsilon=12)
        geometry = [
            mp.Cylinder(radius=R1, center=[0, 0, 0], axis=[0, 0, 1], height=1, material=mat),
            mp.Cylinder(radius=R1, center=[0, 0, 0], axis=[0, 1, 0], height=1, material=mat),
            mp.Cylinder(radius=R1, center=[0, 0, 0], axis=[1, 0, 0], height=1, material=mat),
        ]

        # solve the system
        ms = mpb.ModeSolver(
            num_bands=8,
            geometry_lattice=mp.Lattice(basis1=[1, 0, 0], basis2=[0, 1, 0], basis3=[0, 0, 1],
                size=[1, 1, 1]),
            geometry=geometry,
            resolution=16,
        )
        ms.init_params(p=mp.ALL, reset_fields=true)

        # obtain the symmetry vectors of the bands computed above
        sg_num = 221
        symvecs, topologies = obtain_symmetry_vectors(ms, sg_num)

        for m in symvecs

            for μᴸ in 1:8
                brs = calc_bandreps(sg_num)
                candidatesv = find_bandrep_decompositions(m, brs, μᴸ_min=μᴸ)

                # we should find at least one decomposition
                isempty(candidatesv) && continue

                for candidates in candidatesv

                    # the decomposition should match the symmetry vector
                    @test !isnothing(candidates.longitudinal)
                    @test !isnothing(candidates.apolarv)

                    for nᵀ⁺ᴸ in candidates.apolarv
                        vᵀ = SymmetryVector(nᵀ⁺ᴸ - candidates.longitudinal) # SymVec of nᵀ

                        @test occupation(m) == occupation(vᵀ)
                        @test irreps(m) == irreps(vᵀ)

                        for (i, mult) in enumerate(multiplicities(m))
                            klabel(irreps(m)[i][1]) == "Γ" && continue

                            @test mult == multiplicities(vᵀ)[i]
                        end
                    end

                    # the decompositions should be physical 
                    @test !isnothing(candidates.ps)
                    for p in candidates.ps
                        @test isinteger.(p)
                    end
                end
            end
        end
    end # SG 221
end # EBR decomposition