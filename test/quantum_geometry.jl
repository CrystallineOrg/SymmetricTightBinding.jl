using Test
using Crystalline
using SymmetricTightBinding
using SymmetricTightBinding: solve # not exported
using LinearAlgebra: det, eigvals, tr, Hermitian

# positive semi-definiteness check (`Matrix`, since StaticArrays' analytic 3×3 eigensolver is
# only ~√eps-accurate for the near-degenerate small eigenvalues occurring here)
ispsd(g) = all(≥(-1e-10 * maximum(abs, g)), eigvals(Hermitian(Matrix(g))))

# Reference implementation of the quantum metric, from central finite differences of the band
# projector P(k) = ∑_{n ∈ ns} |uₙₖ⟩⟨uₙₖ| via the identity gᵢⱼ = ½Re Tr[∂ᵢP ∂ⱼP]. Since `P` is
# gauge invariant, this avoids the eigenvector phase ambiguity that would afflict a naive
# finite-differencing of the eigenvectors; it is also entirely independent of the
# sum-over-states machinery that it is used to check here.
function finite_difference_metric(ptbm, k, ns; δ = 1e-5)
    D = length(k)
    P = k′ -> (u = solve(ptbm, k′; bloch_phase=Val(false))[2][:, ns]; u * u')
    ∂Ps = map(1:D) do i
        δk = [j == i ? δ : 0.0 for j in 1:D]
        (P(k .+ δk) - P(k .- δk)) / 2δ
    end
    return [real(tr(∂Ps[i] * ∂Ps[j])) / 2 for i in 1:D, j in 1:D]
end

@testset "Quantum metric" begin
    # a 2D, 2-band model in p2 without TR: has both a nonzero metric and nonzero curvature
    brs₂ = bandreps(2, Val(2); timereversal=false)
    cbr₂ = @composite brs₂[1] + brs₂[3] # (1d|A) + (1c|A)
    ptbm₂ = tb_hamiltonian(cbr₂, [[0,0], [1,0]])([0, .1, 0, 0, .2, 0, 0, .3])
    ks₂ = [[0.1, 0.2], [0.3, -0.4], [-0.25, 0.15]]

    @testset "Symmetry & positive semi-definiteness" begin
        for k in ks₂, n in 1:2
            g = quantummetric(ptbm₂, k, n)
            @test g == transpose(g) # symmetric (exactly, by construction)
            @test ispsd(g)
        end
    end

    @testset "Agreement with finite-difference projector formula (2D)" begin
        for k in ks₂, n in 1:2
            g = quantummetric(ptbm₂, k, n)
            @test g ≈ finite_difference_metric(ptbm₂, k, n) rtol=1e-6
        end

        # also for a model with more than two bands, so that the sum over intermediate
        # states involves more than a single term (a 4-band model in p4mm)
        brs = bandreps(11, Val(2); timereversal=false) # p4mm
        cbr = @composite brs[1] + brs[5] + brs[10] # (2c|A₁) + (1b|A₁) + (1a|A₁), 4 bands
        tbm = tb_hamiltonian(cbr, [[0,0], [1,0]])
        ptbm = tbm(collect(range(0.3, 1.1, length(tbm))))
        for k in ks₂, n in 1:4
            g = quantummetric(ptbm, k, n)
            @test g ≈ finite_difference_metric(ptbm, k, n) rtol=1e-6
            @test ispsd(g)
        end

        # band multiplets: the sum over intermediate states must exclude all of `ns`, which
        # `finite_difference_metric` checks independently via ½Re Tr[∂ᵢP ∂ⱼP] with `P` the
        # multiplet projector
        for k in ks₂, ns in (1:2, 2:3, 3:4, 1:3)
            g = quantummetric(ptbm, k, ns)
            @test g ≈ finite_difference_metric(ptbm, k, ns) rtol=1e-6
            @test ispsd(g)
        end

        # The intra-multiplet terms omitted by the `m ∉ ns` restriction cancel pairwise under
        # `Im` but add under `Re`. So the Berry curvature is exactly additive over the
        # multiplet, while the metric is not - a sharp, integration-free probe of the
        # restriction on both sides.
        for k in ks₂, ns in (1:2, 2:3, 1:3)
            @test berrycurvature(ptbm, k, ns) ≈ sum(berrycurvature(ptbm, k, n) for n in ns)
            @test quantummetric(ptbm, k, ns) ≉ sum(quantummetric(ptbm, k, n) for n in ns)
        end

        # over all bands there are no intermediate states left, so both vanish
        for k in ks₂
            @test quantummetric(ptbm, k, 1:4) ≈ zeros(2, 2) atol=1e-12
            @test berrycurvature(ptbm, k, 1:4) ≈ 0.0 atol=1e-12
        end
    end

    @testset "Agreement with finite-difference projector formula (3D)" begin
        brs = bandreps(2, Val(3); timereversal=false) # P-1 without TR
        cbr = @composite brs[1] + brs[end] # (1d|A) + (1a|B)
        tbm = tb_hamiltonian(cbr, [[1,0,0], [0,1,0], [0,0,1]])
        ptbm = tbm([-0.70, 0.41, -0.89, -0.01, 0.97, 0.41,
                     1.14, -1.22, -0.39, -0.62, 1.03, 1.64])
        for k in ([0.2, -0.7, 0.3], [-0.1, 0.15, 0.42]), n in 1:2
            g = quantummetric(ptbm, k, n)
            @test g == transpose(g)
            @test ispsd(g)
            @test g ≈ finite_difference_metric(ptbm, k, n) rtol=1e-5
        end
    end

    @testset "Consistency with `berrycurvature`" begin
        # positive semi-definiteness of the quantum geometric tensor Q = g - (i/2)Ω implies
        # det g ≥ Ω²/4 in 2D; for a two-band model, Q is rank-1 and the bound is saturated.
        # This ties the (independently validated) metric to the Berry curvature, and would
        # fail under a relative rescaling or a coordinate-convention mismatch of the two
        for k in ks₂, n in 1:2
            g = quantummetric(ptbm₂, k, n)
            Ω = berrycurvature(ptbm₂, k, n)
            @test det(g) ≈ Ω^2/4
        end
    end

    @testset "`quantumgeometrictensor`: Q = g - (i/2)Ω" begin
        # 2D, incl. a multiplet, in the 4-band p4mm model
        brs = bandreps(11, Val(2); timereversal=false) # p4mm
        cbr = @composite brs[1] + brs[5] + brs[10]
        tbm = tb_hamiltonian(cbr, [[0,0], [1,0]])
        ptbm = tbm(collect(range(0.3, 1.1, length(tbm))))
        for k in ks₂, ns in (1, 2, 1:2, 2:3)
            Q = quantumgeometrictensor(ptbm, k, ns)
            @test Q == Q' # Hermitian (exactly, by construction)
            @test ispsd(Q)
            @test real(Q) == quantummetric(ptbm, k, ns)
            @test -2imag(Q[1,2]) ≈ berrycurvature(ptbm, k, ns)
        end

        # 3D
        brs = bandreps(2, Val(3); timereversal=false) # P-1 without TR
        cbr = @composite brs[1] + brs[end] # (1d|A) + (1a|B)
        tbm = tb_hamiltonian(cbr, [[1,0,0], [0,1,0], [0,0,1]])
        ptbm = tbm([-0.70, 0.41, -0.89, -0.01, 0.97, 0.41,
                     1.14, -1.22, -0.39, -0.62, 1.03, 1.64])
        k = [0.2, -0.7, 0.3]
        Q = quantumgeometrictensor(ptbm, k, 1)
        @test Q == Q'
        @test real(Q) == quantummetric(ptbm, k, 1)
        @test -2imag.([Q[2,3], Q[3,1], Q[1,2]]) ≈ berrycurvature(ptbm, k, 1)
    end

    @testset "Nonzero metric with vanishing Berry curvature (PT symmetry)" begin
        # SG 2 (P-1) with TR: both P and T, so Ω ≡ 0 - but the metric survives
        brs = bandreps(2, Val(3); timereversal=true)
        cbr = @composite brs[1] + brs[end] # (1d|A) + (1a|B)
        tbm = tb_hamiltonian(cbr, [[1,0,0], [0,1,0], [0,0,1]])
        ptbm = tbm(collect(range(-0.9, 1.3, length(tbm))))

        k = [0.2, -0.7, 0.3]
        g = quantummetric(ptbm, k, 1)
        @test all(≈(0.0; atol=1e-12), berrycurvature(ptbm, k, 1)) # Ω ≡ 0 ...
        @test maximum(abs, g) > 1e-2                              # ... but g ≠ 0
        @test g ≈ finite_difference_metric(ptbm, k, 1) rtol=1e-5
    end

    @testset "Symmetry" begin
        # `ptbm₂` (p2) retains inversion, so g(k) = g(-k); cf. the analogous test of the
        # evenness of the Berry curvature in `test/berry.jl`
        for k in ks₂
            @test quantummetric(ptbm₂, k, 1) ≈ quantummetric(ptbm₂, -k, 1) atol=1e-10
        end

        # under a point-group operation, g transforms as a rank-2 tensor: g(Rk) = R⁻ᵀ g(k) R⁻¹
        brs = bandreps(11, Val(2); timereversal=false) # p4mm
        cbr = @composite brs[1] + brs[5] # (2c|A₁) + (1b|A₁), 3 bands
        tbm = tb_hamiltonian(cbr, [[0,0], [1,0]])
        ptbm = tbm(collect(range(0.4, 1.2, length(tbm))))
        R = [0 -1; 1 0] # C₄ in reduced coordinates (a square lattice is self-dual)
        for k in ks₂
            g = quantummetric(ptbm, k, 1)
            @test quantummetric(ptbm, R*k, 1) ≈ inv(R)' * g * inv(R) atol=1e-10
        end
    end

    @testset "1D" begin
        brs = bandreps(2, Val(1)) # 1D SG 2 (p-1)
        cbr = @composite brs[1] + brs[end]
        tbm = tb_hamiltonian(cbr, [[0], [1]])
        ptbm = tbm(collect(range(-0.7, 1.1, length(tbm))))
        for k in ([0.13], [-0.32])
            g = quantummetric(ptbm, k, 1)
            @test size(g) == (1, 1)
            @test g ≈ finite_difference_metric(ptbm, k, 1) rtol=1e-6
        end
    end

    @testset "Errors and edge cases" begin
        # a band degeneracy is an error (as for `berrycurvature`)
        brs = bandreps(16, Val(2); timereversal=false) # p6
        cbr = @composite brs[3] # (2b|A), i.e., graphene
        ptbm = tb_hamiltonian(cbr, [[0,0], [1,0]])([0.5, 1.0, 0.0, 0.0, 0.0])
        k_K = [1/3, 1/3] # Dirac point: bands 1 and 2 are degenerate here
        @test_throws ErrorException quantummetric(ptbm, k_K, 1)
        @test_throws ErrorException berrycurvature(ptbm, k_K, 1)
        # ... but the multiplet spanning the degeneracy is well-defined (and here zero, since
        # the model has no further bands to act as intermediate states)
        @test quantummetric(ptbm, k_K, 1:2) ≈ zeros(2, 2) atol=1e-12
        @test berrycurvature(ptbm, k_K, 1:2) ≈ 0.0 atol=1e-12

        # a degeneracy at E ≈ 0 must also be caught (i.e., not only relative to |E|): the
        # Lieb lattice (p4mm, (2c|A₁) + (1a|A₁)) has a flat band at E = 0 that touches the
        # lower band at M
        brs = bandreps(11, Val(2)) # p4mm
        ptbm_L = tb_hamiltonian(@composite(brs[1] + brs[10]), [[0,0]])([0, 0, 1, 1])
        k_M = [1/2, 1/2]
        @test_throws ErrorException quantummetric(ptbm_L, k_M, 2)
        @test_throws ErrorException berrycurvature(ptbm_L, k_M, 2)
        # ... while the multiplet 1:2 is well-defined, and nonzero, with band 3 outside of it
        g = quantummetric(ptbm_L, k_M, 1:2)
        @test g ≈ finite_difference_metric(ptbm_L, k_M, 1:2) rtol=1e-6
        @test maximum(abs, g) > 1

        # the quantum geometry is defined only for Hermitian models
        ptbm_NH = tb_hamiltonian(cbr₂, [[0,0], [1,0]], Val(NONHERMITIAN))(rand(16))
        @test_throws ErrorException quantummetric(ptbm_NH, [0.1, 0.2], 1)
        @test_throws ErrorException berrycurvature(ptbm_NH, [0.1, 0.2], 1)

        @test_throws ErrorException quantummetric(ptbm₂, [0.1, 0.2], 3)   # out of range
        @test_throws ErrorException quantummetric(ptbm₂, [0.1, 0.2], 0)   # out of range
        @test_throws ErrorException quantummetric(ptbm₂, [0.1, 0.2], 1:3) # out of range
        @test_throws ErrorException quantummetric(ptbm₂, [0.1, 0.2], 2:1) # empty
    end
end
