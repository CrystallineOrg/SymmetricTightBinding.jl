using Test
using SymmetricTightBinding

@testset "Berry curvature" begin
    @testset "Berry phase of π around a perturbatively gapped Dirac point" begin
        brs = calc_bandreps(16, Val(2); timereversal=false)
        # pick the (2b|A) EBR, which is just the usual graphene EBR - but now since we break
        # TR, the EBR is not intrinsically connected
        cbr = @composite brs[3] # (2b|A)
        tbm = tb_hamiltonian(cbr, [[0,0],[1,0]])
        # 2nd term of `tbm` is usual graphene nearest-neighbor hopping; 4th term is Haldane-
        # like TR-breaking from next-nearest neighbors (other terms irrelevant here)
        tbm = tbm[[2, 4]]
        ptbm = tbm([1.0, 0.001]) # very small TR-breaking mass term, opening tiny gap at K
                                 # with very localized Berry curvature peak at K
        
        # now sum the Berry curvature over a small square around K-point: should be ~π
        # (it won't be exactly π, since the gap is finite, our integration domain is finite,
        # our sampling is finite, etc.: but should be close). Note that K is at [1/3, 1/3]
        w = 0.05
        δ = 4e-4
        ks = range(1/3-w, 1/3+w; step=δ)[1:end-1] # 250 k-points
        ϕ = sum(berrycurvature(ptbm, [k1, k2], 1) for k1 in ks, k2 in ks; init=0.0) * δ^2
        @test isapprox(ϕ, π; atol=1e-1)

        # test also that Berry curvature is even under inversion (k → -k), since model
        # retains inversion (2-fold rotation) symmetry
        k_samples = [[0.1, 0.2], [0.3, -0.4], [-0.25, 0.15], [1/3, 1/3], [0.0, 0.0]]
        for k in k_samples
            C₊ = berrycurvature(ptbm, k, 1)
            C₋ = berrycurvature(ptbm, -k, 1)
            @test isapprox(C₊, C₋; atol=1e-10)
        end
    end

    @testset "Berry curvature is odd in under TR symmetry with P breaking" begin
        brs = calc_bandreps(13, Val(2); timereversal=true) # p3
        cbr = @composite brs[3] + brs[1] # (1c|A) + (1b|A)
        tbm = tb_hamiltonian(cbr, [[0,1]])
        # terms 1 and 2 break P; term 3 is graphene **next-**nearest-neighbor hopping
        ptbm = tbm([-.2, .2, .5])

        k_samples = [[0.1, 0.2], [0.3, -0.4], [-0.25, 0.15], [1/3, 1/3], [0.0, 0.0]]
        for k in k_samples
            C₊ = berrycurvature(ptbm, k, 1)
            C₋ = berrycurvature(ptbm, -k, 1)
            @test isapprox(C₊, -C₋; atol=1e-10)
            if iszero(k)
                @test isapprox(C₊, 0.0; atol=1e-10)
            else
                @test abs(C₊) > 1e-2
            end
        end
    end

    @testset "Berry curvature in 3D" begin
        # example in which the Berry curvature must be completely zero: both TR + P
        brs = calc_bandreps(2, Val(3); timereversal = true) # P-1 with TR
        cbr = @composite brs[1] + brs[end] # (1d|A) + (1a|B)
        tbm = tb_hamiltonian(cbr, [[1,0,0], [0,1,0], [0,0,1]])

        ptbm = tbm(randn(length(tbm)))
        k = [0.2, -0.7, 0.3]
        Ω = berrycurvature(ptbm, k, 1)
        @test all(≈(0.0), Ω)

        # example where the berry curvature should not vanish: only P
        brs = calc_bandreps(2, Val(3); timereversal = false) # P-1 without TR
        cbr = @composite brs[1] + brs[end] # (1d|A) + (1a|B)
        tbm = tb_hamiltonian(cbr, [[1,0,0], [0,1,0], [0,0,1]])

        ptbm = tbm([-0.70, 0.41, -0.89, -0.01, 0.97, 0.41, 1.14, -1.22, -0.39, -0.62, 1.03, 1.64])
        Ω = berrycurvature(ptbm, k, 1)
        @test Ω ≈ [-6.60267522691563, 6.064530814892586, 1.5689887262383664]

        # example with nonzero Chern number on a k₃ plane: C₃(k₃ = 0.1) = -1
        ptbm = tbm([-0.70, 0.41, -0.89, -0.01, 0.97, 0.41, 1.14, -1.22, -0.39, -0.62, 1.03, 1.64])
        k3 = 0.1
        Nk = 101
        k12s = range(-0.5, 0.5, Nk+1)[2:end]
        Φ = sum(berrycurvature(ptbm, [k1, k2, k3], 1)[3] for k1 in k12s, k2 in k12s; init=0.0)
        C₃ = Φ / (Nk^2 * 2π)
        @test isapprox(C₃, -1.0; atol=1e-2)
    end
end

@testset "Chern numbers" begin
    @testset "Haldane model (p3 w/o TR)" begin
        sgnum = 13 # p3

        brs = calc_bandreps(sgnum, Val(2); timereversal=false)
        # swap (1b|A) and (1c|A) so 1b falls in first orbital (so that Haldane's A position
        # is our "first orbital", and B position our "second orbital"; just makes comparison
        # easier)
        brs[4], brs[1] = brs[1], brs[4]

        # in graphene with a mass staggering term, the relevant EBR is (1b|A) + (1c|A),
        # since the 2b position in plane group p6(mm) subduces to 1b and 1c in p3
        cbr = @composite brs[1] + brs[4] # (1b|A) + (1c|A)
        tbm = tb_hamiltonian(cbr, [[0,0], [0,1]])

        # the below creates a Haldane model with mass term `m`, nearest-neighbor "graphene"
        # hopping `t₁`, and next-nearest-neighbor complex hopping with Haldane phase pattern
        # `t₂ exp(±iϕ)`
        haldane_model(t₁, m, t₂, ϕ) = tbm([m, t₂*cos(ϕ), t₂*sin(ϕ), -m, t₂*cos(ϕ), -t₂*sin(ϕ), t₁, 0])

        # below is a manual version of the Haldane model Hamiltonian, for testing that 
        # `haldane_model` indeed is the correct parameterization (taken from my notes on the
        # Haldane model, Eq. (4))
        function manual_haldane_model(k, t₁, m, t₂, ϕ)
            # constant angles / vectors
            θ₁ = (2*π/3)*1 - π/6; θ₂ = (2*π/3)*2 - π/6; θ₃ = (2*π/3)*3 - π/6
            ξ₁ = (2π/3)*(1-1); ξ₂ = (2π/3)*(2-1); ξ₃ = (2π/3)*(3-1)
            a₁ = -[cos(θ₁), sin(θ₁)]/√(3); a₂ = -[cos(θ₂), sin(θ₂)]/√(3); a₃ = -[cos(θ₃), sin(θ₃)]/√(3)
            b₁ = -[cos(ξ₁), sin(ξ₁)]; b₂ = -[cos(ξ₂), sin(ξ₂)]; b₃ = -[cos(ξ₃), sin(ξ₃)]

            # convert input k-point to cartesian coords
            kc = cartesianize(k, reciprocalbasis(directbasis(13, Val(2))))

            H = (
                t₁ * (cos(dot(kc, a₁)) + cos(dot(kc, a₂)) + cos(dot(kc, a₃))) .* [0 1; 1 0] +
                t₁ * (sin(dot(kc, a₁)) + sin(dot(kc, a₂)) + sin(dot(kc, a₃))) .* [0 -1im; 1im 0] +
                2 * t₂* cos(ϕ) * (cos(dot(kc, b₁)) + cos(dot(kc, b₂)) + cos(dot(kc, b₃))) .* [1 0; 0 1] +
                2 * t₂* sin(ϕ) * (sin(dot(kc, b₁)) + sin(dot(kc, b₂)) + sin(dot(kc, b₃))) .* [1 0; 0 -1] +
                m .* [1 0; 0 -1]
            )
            return H
        end

        # test that our `haldane_model` function indeed gives Haldane's matrix expression
        k_samples = [rand(2) for _ in 1:10] # random k-points for testing
        t₁, m, t₂, ϕ = 0.721, -0.5213, 0.21231, π/3
        for k in k_samples
            H_ptbm = haldane_model(t₁, m, t₂, ϕ)(k)
            H_original = manual_haldane_model(k, t₁, m, t₂, ϕ)

            @test isapprox(H_ptbm, H_original)
        end

        # test that we also get the right Chern numbers
        Nk = 31

        ptbm1 = haldane_model(1.0, 0.5, 0.3, π/2) # must give C = 1 (-1) for 1st (2nd) band
        @test chern_fukui(ptbm1, 1, Nk) == 1
        @test chern(ptbm1, 1, Nk) ≈ 1 atol=3e-2
        @test chern_fukui(ptbm1, 2, Nk) == -1
        @test chern(ptbm1, 2, Nk) ≈ -1 atol=3e-2

        ptbm2 = haldane_model(1.0, 1.5, 0.3, -π/2) # must give C = -1 (+1) for 1st (2nd) band
        @test chern_fukui(ptbm2, 1, Nk) == -1
        @test chern(ptbm2, 1, Nk) ≈ -1 atol=3e-2
        @test chern_fukui(ptbm2, 2, Nk) == 1
        @test chern(ptbm2, 2, Nk) ≈ 1 atol=3e-2

        ptbm3 = haldane_model(1.0, 0.9*3√3, 0.3, 0.1) # must give C = 0 for both bands
        @test chern_fukui(ptbm3, 1, Nk) == 0
        @test chern(ptbm3, 1, Nk) ≈ 0 atol=3e-2
        @test chern_fukui(ptbm3, 2, Nk) == 0
        @test chern(ptbm3, 2, Nk) ≈ 0 atol=3e-2

        # test non-Abelian Fukui method
        @test chern_fukui(ptbm1, [1:1, 2:2], Nk) == [1, -1]
        @test chern_fukui(ptbm1, [1:2], Nk) == [0]
        @test chern_fukui(ptbm1, [1, 2], Nk) == chern_fukui(ptbm1, 1:2, Nk) == only(chern_fukui(ptbm1, [1:2], Nk)) == 0
        @test chern_fukui(ptbm2, [1], Nk) == chern_fukui(ptbm2, 1, Nk)
    end

    @testset "An example in p4" begin
        brs = calc_bandreps(10, Val(2); timereversal=false);
        cbr = @composite brs[9] + brs[10] #  (1a|²E) + (1a|¹E) (2 bands)
        tbm = tb_hamiltonian(cbr, [[0,0], [1,0], [1,1]])

        # construct a model with Chern number -2 (since it has interior Dirac points for
        # zero time-reversal breaking)
        ptbm = tbm([0,1,0,0,0,0,1,0,1,.1]) # 10th term breaks TR (8th also, but set to 0)

        @test chern_fukui(ptbm, 1, 31) == -2
    end
end

