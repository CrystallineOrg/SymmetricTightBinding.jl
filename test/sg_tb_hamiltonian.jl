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
