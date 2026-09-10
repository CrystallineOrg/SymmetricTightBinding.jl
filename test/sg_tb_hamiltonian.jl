using Test
using SymmetricTightBinding
using Crystalline
using LinearAlgebra

# reproducible, generic coefficients
_sg_coefficients(n) = [0.3*cospi(0.73*k) for k in 1:n]

@testset "TB examples in space groups" begin
    @testset "SG 2 (P-1), 3D, single-site EBR" begin
        brs = calc_bandreps(2, Val(3))
        cbr = @composite brs[1]
        tbm = tb_hamiltonian(cbr, [[0, 0, 0], [1, 0, 0]])
        @test length(tbm) > 0

        ptbm = tbm(_sg_coefficients(length(tbm)))
        for k in [[0.1, 0.2, 0.3], [0.4, 0.15, 0.35], [0.0, 0.0, 0.0]]
            H = ptbm(k)
            @test H ≈ H'  # Hermiticity
        end
    end

    @testset "SG 16 (P222), 3D" begin
        brs = calc_bandreps(16, Val(3))
        cbr = @composite brs[1]
        tbm = tb_hamiltonian(cbr, [[0, 0, 0]])

        ptbm = tbm(_sg_coefficients(length(tbm)))
        # H(k) = H*(-k) from TRS
        k = [0.1, 0.2, 0.3]
        H_k = copy(ptbm(k))
        H_mk = copy(ptbm(-k))
        @test H_k ≈ conj.(H_mk)
    end

    @testset "SG 225 (Fm-3m), 3D" begin
        brs = calc_bandreps(225, Val(3))
        cbr = @composite brs[1]
        tbm = tb_hamiltonian(cbr, [[0, 0, 0]])
        @test tbm.N > 0

        if length(tbm) > 0
            ptbm = tbm(_sg_coefficients(length(tbm)))
            H_Γ = ptbm([0.0, 0.0, 0.0])
            @test H_Γ ≈ H_Γ'
        end
    end

    @testset "1D: SG 2 (p-1)" begin
        brs = calc_bandreps(2, Val(1))
        cbr = @composite brs[1]
        tbm = tb_hamiltonian(cbr, [[0], [1]])

        ptbm = tbm(_sg_coefficients(length(tbm)))
        # 1D model: spectrum should be periodic in k
        es_0 = spectrum(ptbm, [0.0])
        es_1 = spectrum(ptbm, [1.0])
        @test es_0 ≈ es_1  atol=1e-10
    end

    @testset "Multi-EBR composite, SG 47 (Pmmm)" begin
        # SG 47: all Wyckoff positions are special (no free parameters)
        brs = calc_bandreps(47, Val(3))
        if length(brs) ≥ 2
            cbr = @composite brs[1] + brs[2]
            tbm = tb_hamiltonian(cbr, [[0, 0, 0]])
            @test tbm.N == occupation(brs[1]) + occupation(brs[2])
        end
    end
end

# ---------------------------------------------------------------------------------------- #

@testset "Diagonal-block terms are not double-counted (issue #131)" begin
    # `evaluate_tight_binding_term!` adds the hermiticity-related block only for *off*-
    # diagonal blocks: a diagonal block is already covered in full by the row/column loop.
    # Adding it there too doubled every off-diagonal element of the block while leaving its
    # diagonal alone - a uniform ×2 for terms confined to one triangle, but a *relative*
    # distortion for terms straddling the block diagonal (multi-dimensional site irreps),
    # which then broke the space-group symmetry outright.

    @testset "graphene nearest-neighbor amplitude" begin
        # textbook: with nearest-neighbor hopping t, |H₁₂(k=0)| = 3t and the bands span ±3t
        brs = calc_bandreps(17, Val(2))
        ptbm = tb_hamiltonian((@composite brs[5]), [[0,0]])([0.0, 1.0])
        @test abs(ptbm([0.0, 0.0])[1, 2]) ≈ 3
        Es = reduce(vcat, [spectrum_single_k(ptbm, [k1, k2])
                           for k1 in range(-0.5, 0.5, 25), k2 in range(-0.5, 0.5, 25)])
        @test maximum(Es) ≈ 3 rtol=1e-2
    end

    @testset "H(gk) is isospectral with H(k) at generic k" begin
        # the irrep-based checks of `test/symmetry_analysis.jl` do not catch this, since
        # they only probe high-symmetry k-points: a generic k is needed
        for (sgnum, Dᵛ, idx, op, Rs) in (
                (17,  Val(2), 5,  S"-y,x-y",   [[0,0], [1,0]]),      # (2b|A₁), 1D site irrep
                (11,  Val(2), 1,  S"-y,x",     [[0,0], [1,0]]),      # (2c|A₁)
                (147, Val(3), 14, S"-y,x-y,z", [[0,0,0], [1,0,0]]),  # (1a|Eᵤ), 2D site irrep
                (147, Val(3), 10, S"-y,x-y,z", [[0,0,0], [1,0,0]]))  # (1b|Eᵤ)
            brs = calc_bandreps(sgnum, Dᵛ)
            cbr = CompositeBandRep([n == idx ? 1 : 0 for n in eachindex(brs)], brs)
            tbm = tb_hamiltonian(cbr, Rs)
            ptbm = tbm(_sg_coefficients(length(tbm)))
            k = Dᵛ === Val(3) ? [0.13, 0.27, 0.19] : [0.13, 0.27]
            gk = rotation(op)' \ k # `k` mapped by `op`, i.e. (R⁻¹)ᵀk
            @test sort(spectrum_single_k(ptbm, k)) ≈ sort(spectrum_single_k(ptbm, gk))
        end
    end
end
