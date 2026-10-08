using Test
using SymmetricTightBinding
using SymmetricTightBinding: site_induced_sgrep_excl_phase, site_induced_timereversal_unitary
using Crystalline
using Crystalline: free
using LinearAlgebra

# reproducible, generic coefficients
_spinful_coefficients(n) = [0.3*cospi(0.73*k) for k in 1:n]

function _spinful_cbr(sgnum, idxs; timereversal::Bool = false)
    brs = bandreps(sgnum; spinful = Val(true), timereversal)
    for i in unique(idxs) # pin free parameters, if any
        iszero(free(position(brs[i]))) || pin_free!(brs, i => [0.13, 0.21, 0.17])
    end
    return CompositeBandRep([count(==(n), idxs) for n in eachindex(brs)], brs)
end

@testset "Spinful models without time-reversal" begin
    # (sgnum, EBR indices, Rs): composites test off-diagonal blocks, whose constraints
    # depend on the SU(2) elements of the generators being chosen identically across blocks
    Rs0 = [[0, 0, 0], [1, 0, 0]]
    Rs1 = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
    cases = [ # sgnum, EBR indices, Rs input
        (16,  [1, 2],   Rs1), # P222: off-diagonal block
        (26,  [1, 3],   Rs1), # Pmc2₁: complex site irreps, but conjugates induce same EBR
        (47,  [1, 2],   Rs1), # Pmmm
        (75,  [3, 7],   Rs1), # P4: free parameters (pinned), complex site irreps
        (143, [1, 4],   Rs1), # P3: real site irreps, complex little group irreps at K & KA
        (166, [13, 19], Rs0), # R-3m: trigonal frame, R centring
        (191, [11, 17], Rs0), # P6/mmm: hexagonal frame
        (225, [12, 16], Rs0), # Fm-3m: F centring, 2D + 4D site irreps
        (227, [13],     Rs0), # Fd-3m: non-symmorphic, F centring
    ]
    ks = [ReciprocalPoint(0.13, 0.27, 0.19), ReciprocalPoint(0.31, -0.07, 0.44)]
    for (sgnum, idxs, Rs) in cases
        cbr = _spinful_cbr(sgnum, idxs)
        @testset "SG $sgnum: $cbr" begin
            ops = primitivize(spacegroup(sgnum, Val(3); spinful = Val(true)))

            # the induced representation is a representation of G/T
            ρs = Dict(g => Matrix(site_induced_sgrep_excl_phase(cbr, g)) for g in ops)
            @test all(Iterators.product(ops, ops)) do (g₁, g₂)
                g₁₂ = ops[findfirst(g -> isapprox(g, g₁ * g₂, nothing, true), ops)]
                ρs[g₁] * ρs[g₂] ≈ ρs[g₁₂]
            end

            tbm = tb_hamiltonian(cbr, Rs)
            @test length(tbm) > 0
            ptbm = tbm(_spinful_coefficients(length(tbm)))
            for k in ks
                Hk = copy(ptbm(k))
                @test Hk ≈ Hk'
                # H(gk) = D_k(g) H(k) D_k(g)† for every operation, not just the generators
                # that the model was built from
                for g in ops
                    D = site_induced_sgrep(ptbm, g)(k)
                    gk = g * k
                    @test ptbm(gk) ≈ D * Hk * D' atol = 1e-12
                end
            end

            # the irreps of the bands reproduce the input EBRs: this checks the SU(2)
            # elements used in the construction against Crystalline's double group irreps
            # NB: skipped for SG 75, whose complex site irreps `symmetry_eigenvalues`
            #     conjugates (see `[⚠️ phase]` and issue #137), identifying the bands as
            #     the EBRs induced from the conjugate site irreps (also in spinless models
            #     without TR)
            @test sum(collect_compatible(ptbm)) == SymmetryVector(cbr) skip = sgnum == 75
        end
    end
end

@testset "Spinful models with time-reversal" begin
    Rs0 = [[0, 0, 0], [1, 0, 0]]
    Rs1 = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
    cases = [ # sgnum, EBR indices, Rs input
        (2,   [1],    Rs1), # P-1: (1h|AᵤˢAᵤˢ), inversion ⇒ Kramers degeneracy at all k
        (16,  [1, 2], Rs1), # P222: off-diagonal block, no inversion
        (75,  [1, 2], Rs1), # P4: glued complex site irreps (e.g., ¹E₁ˢ²E₁ˢ)
        (143, [1, 2], Rs1), # P3: trigonal frame, (1c|EˢEˢ) + (1c|¹Eˢ²Eˢ)
        (166, [1, 2], Rs0), # R-3m: R centring
        (191, [1],    Rs0), # P6/mmm: hexagonal frame
        (225, [1],    Rs0), # Fm-3m: F centring
        (227, [1],    Rs0), # Fd-3m: non-symmorphic, F centring
    ]
    ks = [ReciprocalPoint(0.13, 0.27, 0.19), ReciprocalPoint(0.31, -0.07, 0.44)]
    for (sgnum, idxs, Rs) in cases
        cbr = _spinful_cbr(sgnum, idxs; timereversal = true)
        @testset "SG $sgnum: $cbr" begin
            ops = primitivize(spacegroup(sgnum, Val(3); spinful = Val(true)))
            tbm = tb_hamiltonian(cbr, Rs)
            @test length(tbm) > 0
            ptbm = tbm(_spinful_coefficients(length(tbm)))
            Γ = site_induced_timereversal_unitary(cbr)
            has_inversion = any(g -> rotation(g) == -I, ops)
            for k in ks
                Hk = copy(ptbm(k))
                @test Hk ≈ Hk'
                # H(-k) = Γ H*(k) Γ†, with `Γ = 𝟙 ⊗ J` (`J = iσʸ ⊗ 𝟙ₙ`) on each site
                @test ptbm(-k) ≈ Γ * conj(Hk) * Γ' atol = 1e-12
                for g in ops
                    D = site_induced_sgrep(ptbm, g)(k)
                    @test ptbm(g * k) ≈ D * Hk * D' atol = 1e-12
                end
                # Kramers degeneracy at generic k iff inversion is present (since 𝒯² = -1)
                E = eigvals(Hermitian(Hk))
                @test all(i -> isapprox(E[2i-1], E[2i]; atol = 1e-9), 1:length(E)÷2) == has_inversion
            end
            @test sum(collect_compatible(ptbm)) == SymmetryVector(cbr)
        end
    end

    # hand-check: in P-1, a single Kramers pair (1h|AᵤˢAᵤˢ) with inversion and time-reversal
    # has, per hopping vector, a single real hopping amplitude (∝ 𝟙): i.e., 1 on-site term
    # and 1 term for each of the 3 nearest-neighbor hopping vectors
    @test length(tb_hamiltonian(_spinful_cbr(2, [1]; timereversal = true), Rs1)) == 4
end

@testset "`isspinful` for tight-binding models" begin
    Rs = [[0, 0, 0]]
    for spinful in (Val(false), Val(true))
        brs = bandreps(16; spinful, timereversal = false)
        cbr = @composite brs[1] + brs[2]
        tbm = tb_hamiltonian(cbr, Rs)
        ptbm = tbm(_spinful_coefficients(length(tbm)))
        ctbm = tbm + tb_hamiltonian(cbr, Rs, Val(ANTIHERMITIAN))
        pctbm = ctbm(_spinful_coefficients(length(ctbm.tbm_h)),
                     _spinful_coefficients(length(ctbm.tbm_a)))
        spinful_value = spinful === Val(true) ? true : false
        for x in (tbm[1].block, tbm[1], tbm, ptbm, ctbm, pctbm)
            @test isspinful(x) == spinful_value
            @test isspinful(typeof(x)) == spinful_value
        end
    end
end

@testset "Spinful models: unsupported input" begin
    cbr = _spinful_cbr(16, [1, 2])
    tbm = tb_hamiltonian(cbr, [[0, 0, 0]])
    @test_throws "requires a double group operation" site_induced_sgrep(tbm, S"-x,-y,z")
end

