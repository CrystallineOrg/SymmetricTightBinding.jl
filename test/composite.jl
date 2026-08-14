using Test
using SymmetricTightBinding
using Crystalline
using LinearAlgebra: ishermitian
using SymmetricTightBinding: solve, TightBindingTerm # (`orbital_count` &c are exported)

# `test_show` & `test_tp_show`; guarded so this file also runs standalone
isdefined(@__MODULE__, :test_show) || include("test_utils.jl")

@testset "CompositeTightBindingModel" begin

    # 1D non-Hermitian SSH model, (1b|A′) ⊕ (1a|A′) in 1D SG 2, as in test/nonhermitian.jl
    brs = calc_bandreps(2, Val(1))
    cbr = @composite brs[1] + brs[3]
    Rs = [[0,], [1,]]
    tbm_h = tb_hamiltonian(cbr, Rs, Val(HERMITIAN))
    tbm_a = tb_hamiltonian(cbr, Rs, Val(ANTIHERMITIAN))
    ctbm = tbm_h + tbm_a
    Nʰ, Nᵃ = length(tbm_h), length(tbm_a)
    @test (Nʰ, Nᵃ) == (5, 1)

    @testset "Construction" begin
        # all spellings are equivalent, and the argument order is immaterial
        for c in (CompositeTightBindingModel(tbm_h, tbm_a),
                  CompositeTightBindingModel(tbm_a, tbm_h),
                  CompositeTightBindingModel{1}(tbm_h, tbm_a),
                  CompositeTightBindingModel{1}(tbm_a, tbm_h),
                  tbm_h + tbm_a,
                  tbm_a + tbm_h)
            @test c.tbm_h === tbm_h && c.tbm_a === tbm_a
        end

        # the two parts must be built on the same composite band representation
        cbr′ = @composite brs[1] + brs[1]
        tbm_a′ = tb_hamiltonian(cbr′, Rs, Val(ANTIHERMITIAN))
        @test_throws ErrorException CompositeTightBindingModel(tbm_h, tbm_a′)
    end

    @testset "Free-parameter count matches the NONHERMITIAN model" begin
        # a composite model is a change of basis of the corresponding non-Hermitian model,
        # so the Hermitian and anti-Hermitian terms must together span the same space
        for (sgnum, Dᵛ, Rs′) in ((2, Val(1), [[0,], [1,]]),
                                 (10, Val(2), [[0,0], [1,0]]),
                                 (13, Val(2), [[0,0], [1,0]]))
            brs′ = calc_bandreps(sgnum, Dᵛ)
            cbr′ = sgnum == 2 ? (@composite brs′[1] + brs′[3]) : (@composite brs′[1])
            n_h = length(tb_hamiltonian(cbr′, Rs′, Val(HERMITIAN)))
            n_a = length(tb_hamiltonian(cbr′, Rs′, Val(ANTIHERMITIAN)))
            n_nh = length(tb_hamiltonian(cbr′, Rs′, Val(NONHERMITIAN)))
            @test n_h + n_a == n_nh
        end
    end

    @testset "AbstractVector interface" begin
        @test length(ctbm) == Nʰ + Nᵃ
        @test size(ctbm) == (Nʰ + Nᵃ,)
        @test eltype(ctbm) == Union{TightBindingTerm{1, HERMITIAN},
                                    TightBindingTerm{1, ANTIHERMITIAN}}
        # Hermitian terms come first, then anti-Hermitian ones
        @test all(i -> ctbm[i] === tbm_h[i], 1:Nʰ)
        @test all(i -> ctbm[Nʰ+i] === tbm_a[i], 1:Nᵃ)
        @test collect(ctbm) == vcat(collect(tbm_h), collect(tbm_a))
        @test_throws BoundsError ctbm[0]
        @test_throws BoundsError ctbm[Nʰ+Nᵃ+1]
    end

    @testset "Indexing with index vectors" begin
        # ranges, gappy index vectors, and masks all give a composite model back, with the
        # requested terms in the requested order
        for idxs in (1:3, 4:6, [1, 3, 6], setdiff(1:6, 3), [6], Int[], 1:6)
            sub = ctbm[idxs]
            @test length(sub) == length(idxs)
            @test all(j -> sub[j] === ctbm[idxs[j]], eachindex(idxs))
            @test length(sub.tbm_h) == count(≤(Nʰ), idxs)
            @test length(sub.tbm_a) == count(>(Nʰ), idxs)
        end

        # logical indexing (note `Bool <: Integer`, so this needs its own method)
        mask = [true, false, false, false, false, true]
        @test collect(ctbm[mask]) == [ctbm[1], ctbm[6]]
        @test length(ctbm[trues(6)]) == 6
        @test length(ctbm[falses(6)]) == 0
        @test_throws BoundsError ctbm[falses(3)]

        # a composite necessarily stores Hermitian terms before anti-Hermitian ones, so a
        # permuted (or repeated) index vector cannot be honoured: it must error rather than
        # silently return a differently-ordered model
        @test_throws ErrorException ctbm[[6, 1]]
        @test_throws ErrorException ctbm[[3, 1]]
        @test_throws ErrorException ctbm[[1, 1]]
        @test_throws BoundsError ctbm[[0, 1]]
        @test_throws BoundsError ctbm[[1, 7]]
    end

    @testset "Accessors" begin
        @test hermiticity(ctbm) == NONHERMITIAN
        @test orbital_count(ctbm) == orbital_count(tbm_h) == 2
        @test orbital_positions(ctbm) == orbital_positions(tbm_h)
        @test CompositeBandRep(ctbm) == cbr
        @test Crystalline.dim(ctbm) == 1
    end

    @testset "Evaluation & Hermitian/anti-Hermitian splitting" begin
        cs_h = [0.3*cospi(0.73*k) for k in 1:Nʰ]
        cs_a = [0.3*cospi(0.73*k) for k in (Nʰ+1):(Nʰ+Nᵃ)]
        pctbm = ctbm(cs_h, cs_a)
        @test hermiticity(pctbm) == NONHERMITIAN
        @test orbital_count(pctbm) == 2
        @test CompositeBandRep(pctbm) == cbr

        # passing the coefficients jointly or separately must agree
        @test pctbm.cs == vcat(cs_h, cs_a)
        @test ctbm(vcat(cs_h, cs_a)).cs == pctbm.cs
        @test_throws ErrorException ctbm(cs_h, [cs_a; 0.0])  # wrong count
        @test_throws ErrorException ctbm(cs_h)               # too few coefficients

        # the split is canonical: at every k, the Hermitian part of H(k) is (H+H′)/2 and
        # the anti-Hermitian part is (H-H′)/2 (NB: `pctbm(k)` aliases scratch, so `copy`)
        ptbm_h, ptbm_a = tbm_h(cs_h), tbm_a(cs_a)
        for k in ([0.0], [0.2], [0.37], [0.5])
            H = copy(pctbm(k))
            Hʰ, Hᵃ = Matrix(ptbm_h(k)), copy(ptbm_a(k))
            @test Hʰ + Hᵃ ≈ H
            @test Hʰ ≈ (H + H') / 2
            @test Hᵃ ≈ (H - H') / 2
        end

        # without anti-Hermitian terms, the model is Hermitian at every k
        pctbm_h = ctbm(cs_h, zeros(Nᵃ))
        @test all(k -> ishermitian(pctbm_h([k])), (0.0, 0.2, 0.37, 0.5))
        # ... and genuinely non-Hermitian once they are switched on
        @test !ishermitian(pctbm([0.2]))

        # spectra of a composite model are complex, and `solve` agrees with `spectrum`
        Es = spectrum(pctbm, [[0.0], [0.2], [0.5]])
        @test eltype(Es) == ComplexF64
        @test size(Es) == (3, 2)
        by_reim = x -> (real(x), imag(x)) # `eigvals!`/`eigen!` need not order identically
        @test sort(Es[2, :]; by = by_reim) ≈ sort(solve(pctbm, [0.2])[1]; by = by_reim)
        @test length(solve(pctbm, [0.2]; bloch_phase = Val(true))[1]) == 2

        # a term-dropped sub-model still evaluates, and still splits canonically
        sub = ctbm[[1, 2, 6]]
        psub = sub([0.3*cospi(0.73*k) for k in 1:3])
        Hˢ = copy(psub([0.2]))
        @test Matrix(sub.tbm_h(psub.cs[1:2])([0.2])) ≈ (Hˢ + Hˢ') / 2
    end

    @testset "Show methods" begin
        test_show(repr(MIME"text/plain"(), ctbm),
        """
        (5+1)-term 2×2 CompositeTightBindingModel{1} over (1b|A′)⊕(1a|A′):
        ┌─ Hermitian
        1. ⎡ 1  │  0 ⎤
        │  ⎢ ───┼─── ⎥
        │  ⎣ 0  │  0 ⎦
        └─ (1b|A′) self-term
        ┌─
        2. ⎡ 𝕖(δ₁)+𝕖(δ₂)  │  0 ⎤
        │  ⎢ ─────────────┼─── ⎥
        │  ⎣ 0            │  0 ⎦
        └─ (1b|A′) self-term:  δ₁=[-1], δ₂=-δ₁
        ┌─
        3. ⎡ 0  │  0 ⎤
        │  ⎢ ───┼─── ⎥
        │  ⎣ 0  │  1 ⎦
        └─ (1a|A′) self-term
        ┌─
        4. ⎡ 0  │  0           ⎤
        │  ⎢ ───┼───────────── ⎥
        │  ⎣ 0  │  𝕖(δ₁)+𝕖(δ₂) ⎦
        └─ (1a|A′) self-term:  δ₁=[-1], δ₂=-δ₁
        ┌─
        5. ⎡ 0              │  𝕖(δ₁)+𝕖(δ₂) ⎤
        │  ⎢ ───────────────┼───────────── ⎥
        │  ⎣ 𝕖(-δ₁)+𝕖(-δ₂)  │  0           ⎦
        └─ (1b|A′)↔(1a|A′):  δ₁=[1/2], δ₂=-δ₁
        ┌─ Anti-Hermitian
        6. ⎡ 0               │  𝕖(δ₁)+𝕖(δ₂) ⎤
        │  ⎢ ────────────────┼───────────── ⎥
        │  ⎣ -𝕖(-δ₁)-𝕖(-δ₂)  │  0           ⎦
        └─ (1b|A′)↔(1a|A′):  δ₁=[1/2], δ₂=-δ₁""")

        @test sprint(summary, ctbm) ==
            "(5+1)-term 2×2 CompositeTightBindingModel{1} over (1b|A′)⊕(1a|A′)"

        pctbm = ctbm([0.3*cospi(0.73*k) for k in 1:6])
        test_show(repr(MIME"text/plain"(), pctbm),
        """
        (5+1)-term 2×2 ParameterizedCompositeTightBindingModel{1} over (1b|A′)⊕(1a|A′) \
        with amplitudes:
         [-0.19839, -0.0376, 0.24812, -0.29057, 0.1362, 0.11044]""")
    end
end
