using Test
using SymmetricTightBinding
using SymmetricTightBinding: sgrep_induced_by_siteir_excl_phase
using Crystalline
using Crystalline: free
using LinearAlgebra

# reproducible, generic coefficients
_spinful_coefficients(n) = [0.3*cospi(0.73*k) for k in 1:n]

function _spinful_cbr(sgnum, idxs)
    brs = bandreps(sgnum; spinful = Val(true), timereversal = false)
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
            ρs = Dict(g => Matrix(sgrep_induced_by_siteir_excl_phase(cbr, g)) for g in ops)
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
                    D = sgrep_induced_by_siteir(ptbm, g)(k)
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

@testset "Spinful models: unsupported input" begin
    cbr = _spinful_cbr(16, [1, 2])
    tbm = tb_hamiltonian(cbr, [[0, 0, 0]])
    @test_throws "requires a double group operation" sgrep_induced_by_siteir(tbm, S"-x,-y,z")
end

