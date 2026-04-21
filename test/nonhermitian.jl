using Test
using SymmetricTightBinding
using Crystalline
using LinearAlgebra: ishermitian, eigvals

@testset "NONHERMITIAN models" begin

    @testset "NH-SSH model (1D SG 2)" begin
        # (1b|A′) ⊕ (1a|A′) in 1D SG 2 (inversion symmetry); intra- and inter-cell hoppings
        brs = calc_bandreps(2, Val(1))
        cbr = @composite brs[1] + brs[3]
        tbm_NH = tb_hamiltonian(cbr, [[0,], [1,]], Val(NONHERMITIAN))
        tbm_H  = tb_hamiltonian(cbr, [[0,], [1,]], Val(HERMITIAN))

        @test hermiticity(tbm_NH) == NONHERMITIAN
        @test hermiticity(tbm_H)  == HERMITIAN

        # NONHERMITIAN: upper & lower off-diagonal blocks are independent → 6 terms
        # HERMITIAN:    lower off-diagonal block follows from upper by conjugation  → 5 terms
        # (all 4 diagonal self-terms are shared; only off-diagonal count differs: 2 vs 1)
        @test length(tbm_NH) == 6
        @test length(tbm_H)  == 5

        # equal off-diagonal coefficients (c₅ = c₆) → H(k) is Hermitian for all k
        ptbm = tbm_NH([0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
        @test all(k -> ishermitian(ptbm([k])), [0.0, 0.25, 0.5])

        # unequal off-diagonal coefficients → genuinely non-Hermitian
        ptbm = tbm_NH([0.0, 0.0, 0.0, 0.0, 1.0, 0.5])
        @test !ishermitian(ptbm([0.3]))

        # spectrum is complex-valued for a NONHERMITIAN model
        @test eltype(spectrum(ptbm, [[0.0], [0.3], [0.5]])) == ComplexF64
    end

    @testset "2D NONHERMITIAN with p4 symmetry (PG 10)" begin
        # 2-site EBR on p4 lattice (C4, no mirrors); NN + on-site hoppings
        # C4 constrains all 4 NN directions into a single orbit; NONHERMITIAN adds
        # an independent lower-triangular off-diagonal block → more free parameters
        brs = calc_bandreps(10, Val(2))
        cbr = @composite brs[1]
        tbm_NH = tb_hamiltonian(cbr, [[0,0], [1,0]], Val(NONHERMITIAN))
        tbm_H  = tb_hamiltonian(cbr, [[0,0], [1,0]], Val(HERMITIAN))

        @test hermiticity(tbm_NH) == NONHERMITIAN
        @test hermiticity(tbm_H)  == HERMITIAN
        # NONHERMITIAN has more free parameters than HERMITIAN (independent off-diagonal blocks)
        @test length(tbm_NH) > length(tbm_H)
        @test length(tbm_NH) - length(tbm_H) == 2  # two extra blocks (one per off-diagonal orbit)
    end

    @testset "NONHERMITIAN without TR allows complex on-site (1D SG 2)" begin
        # single EBR (1b|A′) in 1D SG 2; with TR the on-site energy is forced real by
        # H(k) = H*(−k); without TR it can be complex (uniform gain/loss)
        brs_TR   = calc_bandreps(2, Val(1))
        brs_noTR = calc_bandreps(2, Val(1); timereversal = false)
        cbr_TR   = @composite brs_TR[1]
        cbr_noTR = @composite brs_noTR[1]
        tbm_TR   = tb_hamiltonian(cbr_TR,   [[0,]], Val(NONHERMITIAN))
        tbm_noTR = tb_hamiltonian(cbr_noTR, [[0,]], Val(NONHERMITIAN))

        @test length(tbm_TR)   == 1  # real on-site only (TR forces real)
        @test length(tbm_noTR) == 2  # real on-site + imaginary on-site

        # the imaginary on-site term gives purely imaginary eigenvalues → non-real spectrum
        ptbm = tbm_noTR([0.0, 1.0])  # c₁ = 0 (real part), c₂ = 1 (imaginary part = gain/loss)
        λs = eigvals(ptbm([0.0]))
        @test !all(isreal, λs)
    end

end
