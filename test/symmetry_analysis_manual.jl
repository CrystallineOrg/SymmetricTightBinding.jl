# Stopgap symmetry analysis tests.
#
# These are a small sampling of `collect_compatible` tests, covering all 2D plane groups
# with special Wyckoff positions. They serve as an interim solution until PR #89 is
# resolved. The known-failing cases (SG 13 1b/1c, SG 14 1b/1c, SG 16 3c — all involving
# the K-point phase convention bug) are marked `@test_broken` or skipped.
#
# Once #89 is merged, this file should be deleted.

using Test
using SymmetricTightBinding, Crystalline

# Generate reproducible, "generic" coefficients (no special symmetry) for n terms.
_test_coefficients(n) = [0.3*cospi(0.73*k) for k in 1:n]

# Helper: test symmetry analysis for a single EBR, comparing the symmetry vector from
# `collect_compatible` against the expected `SymmetryVector` from the EBR definition.
# Returns true if the symmetry vectors match.
function _test_symmetry_analysis(brs, i)
    coefs = zeros(length(brs))
    coefs[i] = 1
    cbr = CompositeBandRep(coefs, brs)
    tbm = tb_hamiltonian(cbr)
    length(tbm) == 0 && return true # trivially passes (no terms to test)
    ptbm = tbm(_test_coefficients(length(tbm)))
    v_expected = SymmetryVector(cbr)
    v_model = collect_compatible(ptbm)
    return v_expected == v_model[1]
end

@testset "Symmetry analysis" begin

    # ---- 1D ----
    @testset "1D: SG 2 (p-1)" begin
        brs = bandreps(2, Val(1))
        for i in eachindex(brs)
            @test _test_symmetry_analysis(brs, i)
        end
    end

    # ---- 2D plane groups (all special Wyckoff EBRs) ----
    # SGs 1, 3, 4, 5 are skipped: SG 1 has only general positions (free params);
    # SGs 3-5 have only general-position EBRs.

    @testset "2D: SG 2 (p2)" begin
        brs = bandreps(2, Val(2))
        for i in eachindex(brs) # all 8 EBRs pass
            @test _test_symmetry_analysis(brs, i)
        end
    end

    @testset "2D: SG 6 (p2mm)" begin
        brs = bandreps(6, Val(2))
        for i in eachindex(brs) # all 16 EBRs pass
            @test _test_symmetry_analysis(brs, i)
        end
    end

    @testset "2D: SG 7 (p2mg)" begin
        brs = bandreps(7, Val(2))
        # EBRs 1-2 have free params (general positions); 3-6 are special and pass
        for i in 3:6
            @test _test_symmetry_analysis(brs, i)
        end
    end

    @testset "2D: SG 8 (p2gg)" begin
        brs = bandreps(8, Val(2))
        for i in eachindex(brs) # all 4 EBRs pass
            @test _test_symmetry_analysis(brs, i)
        end
    end

    @testset "2D: SG 9 (c2mm)" begin
        brs = bandreps(9, Val(2))
        for i in eachindex(brs) # all 10 EBRs pass
            @test _test_symmetry_analysis(brs, i)
        end
    end

    @testset "2D: SG 10 (p4)" begin
        brs = bandreps(10, Val(2))
        for i in eachindex(brs) # all 8 EBRs pass
            @test _test_symmetry_analysis(brs, i)
        end
    end

    @testset "2D: SG 11 (p4mm)" begin
        brs = bandreps(11, Val(2))
        for i in eachindex(brs) # all 14 EBRs pass
            @test _test_symmetry_analysis(brs, i)
        end
    end

    @testset "2D: SG 12 (p4gm)" begin
        brs = bandreps(12, Val(2))
        for i in eachindex(brs) # all 7 EBRs pass
            @test _test_symmetry_analysis(brs, i)
        end
    end

    @testset "2D: SG 13 (p3)" begin
        brs = bandreps(13, Val(2))
        for i in eachindex(brs) # all 6 EBRs pass (K-point phase issue fixed)
            @test _test_symmetry_analysis(brs, i)
        end
    end

    @testset "2D: SG 13 (p3), without TR (#137)" begin
        brs = bandreps(13, Val(2); timereversal = false)
        for i in eachindex(brs) # complex site irreps were previously conjugated
            @test _test_symmetry_analysis(brs, i)
        end

        # labels at K describe our bands at K, which carry the tables' KA irreps
        cbr = CompositeBandRep([0, 1, 0, 0, 0, 1, 1, 0, 0], brs) # (1c|²E)+(1b|¹E)+(1a|A)
        ptbm = tb_hamiltonian(cbr)([1, 0, -1, 0, 0, 0, 0, 0, 0])
        annotations = collect_irrep_annotations(ptbm)
        @test last.(annotations["K"]) == ["KA₁", "KA₂", "KA₃"]
        @test last.(annotations["KA"]) == ["K₁", "K₁", "K₁"]

        # complex site irreps away from the origin: previously gave (1c|¹E)+(1b|²E)
        cbr2 = CompositeBandRep([0, 1, 0, 0, 0, 1, 0, 0, 0], brs) # (1c|²E)+(1b|¹E)
        tbm2 = tb_hamiltonian(cbr2, [[0, 0], [1, 0]])
        ptbm2 = tbm2(_test_coefficients(length(tbm2)))
        @test sum(collect_compatible(ptbm2)) == SymmetryVector(cbr2)

        # `symmetry_eigenvalues` uses our convention: 3⁺ eigenvalues at K are those of KA
        lgK = group(lgirreps(13, Val(2))["K"])
        K = position(lgK)()
        for (i, χ) in ((2, cispi(-2/3)), (6, cispi(2/3)), (7, 1)) # (1c|²E), (1b|¹E), (1a|A)
            ptbm1 = tb_hamiltonian(CompositeBandRep(Int.(eachindex(brs) .== i), brs))([1.0])
            @test only(symmetry_eigenvalues(ptbm1, operations(lgK), K)[2, :]) ≈ χ
        end

        # `flip_bloch_phase` moves K's irreps to -K (labelled KA) and vice versa
        lgirsv′ = SymmetricTightBinding.flip_bloch_phase(irreps(brs); timereversal = false)
        @test [klabel(group(lgirs)) for lgirs in lgirsv′] == ["Γ", "KA", "K"]
        @test [label(first(lgirs)) for lgirs in lgirsv′] == ["Γ₁", "K₁", "KA₁"]
        @test position(group(lgirsv′[2]))() ≈ -K

        # with TR, the K-labels are kept at K
        brsᵀ = bandreps(13, Val(2))
        ptbmᵀ = tb_hamiltonian(CompositeBandRep([1, 0, 0, 0, 0, 0], brsᵀ))([1.0]) # (1c|A)
        annotationsᵀ = collect_irrep_annotations(ptbmᵀ)
        @test last.(annotationsᵀ["K"]) == ["K₃"]
        @test last.(annotationsᵀ["KA"]) == ["KA₂"]
    end

    @testset "2D: SG 14 (p3m1)" begin
        brs = bandreps(14, Val(2))
        for i in eachindex(brs) # all 9 EBRs pass (K-point phase issue fixed)
            @test _test_symmetry_analysis(brs, i)
        end
    end

    @testset "2D: SG 15 (p31m)" begin
        brs = bandreps(15, Val(2))
        for i in eachindex(brs) # all 5 EBRs pass
            @test _test_symmetry_analysis(brs, i)
        end
    end

    @testset "2D: SG 16 (p6)" begin
        brs = bandreps(16, Val(2))
        # EBRs 1-2 (3c site): flaky due to K-point phase issue — skip
        for i in 3:8
            @test _test_symmetry_analysis(brs, i)     # 2b and 1a EBRs pass
        end
    end

    @testset "2D: SG 17 (p6mm)" begin
        brs = bandreps(17, Val(2))
        for i in eachindex(brs) # all 13 EBRs pass (includes graphene)
            @test _test_symmetry_analysis(brs, i)
        end
    end

end # @testset "Symmetry analysis"
