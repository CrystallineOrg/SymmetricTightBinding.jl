using TETB, Test
using Crystalline

@testset "EBR decomposition" begin

    @testset "SG #2" begin
        sgnum, D = 2, 3
        brs = calc_bandreps(sgnum, Val(D))
        lgirsv = irreps(brs)
        ms = SymmetryVector{D}[]

        # do some checks for interesting symmetry vectors
        s1 = "[-Γ₁⁺+3Γ₁⁻,2R₁⁻, 2T₁⁻, 2U₁⁻, V₁⁺+V₁⁻, X₁⁺+X₁⁻,Y₁⁺+Y₁⁻,2Z₁⁺]"
        append!(ms, parse(SymmetryVector, s1, lgirsv))

        s2 = "[-Γ₁⁺+3Γ₁⁻, 2R₁⁻, T₁⁺ + T₁⁻, U₁⁺ + U₁⁻, V₁⁺ + V₁⁻, 2X₁⁻, 2Y₁⁻, 2Z₁⁺]"
        append!(ms, parse(SymmetryVector, s2, lgirsv))

        s3 = "[-Γ₁⁺+3Γ₁⁻, R₁⁺+R₁⁻, T₁⁺ + T₁⁻, 2 U₁⁻, 2 V₁⁻, X₁⁺+X₁⁻, 2Y₁⁻, 2Z₁⁺]"
        append!(ms, parse(SymmetryVector, s3, lgirsv))

        s4 = "[-Γ₁⁺+3Γ₁⁻, R₁⁺+R₁⁻, 2 T₁⁻, U₁⁺ + U₁⁻, 2 V₁⁻, 2X₁⁻, Y₁⁺+Y₁⁻, 2Z₁⁺]"
        append!(ms, parse(SymmetryVector, s4, lgirsv))

        s5 = "[-Γ₁⁺+3Γ₁⁻, 2R₁⁻, T₁⁺ + T₁⁻, 2 U₁⁻, 2 V₁⁺, 2X₁⁻, Y₁⁺+Y₁⁻, Z₁⁺+Z₁⁻]"
        append!(ms, parse(SymmetryVector, s5, lgirsv))

        for m in mv
            for μᴸ in 1:8
                brs = calc_bandreps(sgnum)
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
    end # SG 2

end # EBR decomposition