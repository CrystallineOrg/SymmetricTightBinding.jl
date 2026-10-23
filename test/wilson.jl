using Test
using Crystalline
using Crystalline: constant
using SymmetricTightBinding
using LinearAlgebra: dot

# net winding (in units of 2π) of a sequence of phases, unwrapping jumps across ±π
winding(phases) = sum(rem2pi.(diff(phases), RoundNearest)) / 2π

# whether two phases are equal modulo 2π
phase_isapprox(θ, θ′; atol = 1e-8) = isapprox(cis(θ), cis(θ′); atol)

# the EBR whose symmetry vector matches that of the lowest band of `ptbm` (or `nothing`)
function lowest_band_ebr(ptbm, brs)
    sv = first(collect_compatible(ptbm))
    i = findfirst(br -> SymmetryVector(br) == sv, brs)
    return isnothing(i) ? nothing : brs[i]
end

# minimum direct gap between the lowest two bands of a 2D model, over a `Nk×Nk` mesh
function lowest_gap(ptbm; Nk = 60)
    k12s = range(0, 1, Nk+1)[1:end-1]
    Es = spectrum(ptbm, [[k1, k2] for k1 in k12s for k2 in k12s])
    return minimum(Es[:, 2] .- Es[:, 1])
end

@testset "Wilson loops" begin
    # Haldane model, as in test/berry.jl
    brs = bandreps(13, Val(2); timereversal = false) # p3
    brs[4], brs[1] = brs[1], brs[4]
    cbr = @composite brs[1] + brs[4] # (1b|A) + (1c|A)
    tbm = tb_hamiltonian(cbr, [[0,0], [0,1]])
    haldane_model(t₁, m, t₂, ϕ) = tbm([m, t₂*cos(ϕ), -t₂*sin(ϕ), -m, t₂*cos(ϕ), t₂*sin(ϕ), t₁, 0])

    @testset "Winding of Berry phases equals Chern number (Haldane model)" begin
        Nk = 31
        k12s = range(-0.5, 0.5, Nk)
        for (params, Cs) in (((1.0, 0.5, 0.3, π/2), (1, -1)),
                             ((1.0, 1.5, 0.3, -π/2), (-1, 1)),
                             ((1.0, 0.9*3√3, 0.3, 0.1), (0, 0)))
            ptbm = haldane_model(params...)
            for (n, C) in enumerate(Cs)
                # loops along k₂, advanced along k₁
                γs = [wilsonloop(ptbm, [[k1, k2] for k2 in k12s], n) for k1 in k12s]
                @test winding(γs) ≈ C atol = 1e-6
                @test chern_fukui(ptbm, n, Nk) == C # consistency with existing Chern tools
                # reversing the roles of k₁ and k₂ reverses the orientation, flipping the sign
                γs′ = [wilsonloop(ptbm, [[k2, k1] for k2 in k12s], n) for k1 in k12s]
                @test winding(γs′) ≈ -C atol = 1e-6
            end
        end
    end

    # over *all* bands, the link matrices telescope, leaving only the gauge-closing step at
    # the end of the loop: the Wilson loop then reduces to `exp(iG⋅r)`, over the orbital
    # positions r and the loop's (reciprocal lattice) winding G, regardless of mesh or model
    @testset "Full-band Wilson loop reduces to exp(iG⋅r)" begin
        i1a = findfirst(br -> string(br) == "(1a|A)", brs)
        i1b = findfirst(br -> string(br) == "(1b|A)", brs)
        coefs(is...) = [count(==(i), is) for i in eachindex(brs)]
        # orbitals at 1b, and at 1b & 1a: positions chosen so that G⋅r ≢ -G⋅r (mod 1)
        for cbr′ in (CompositeBandRep(coefs(i1b), brs), CompositeBandRep(coefs(i1a, i1b), brs))
            tbm′ = tb_hamiltonian(cbr′, [[0,0], [1,0]])
            ptbm′ = tbm′([0.3*cospi(0.73*k) for k in 1:length(tbm′)])
            N = orbital_count(ptbm′)
            rs = orbital_positions(ptbm′)
            for G in ([0,0], [0,1], [1,0], [1,1], [2,-1])
                ks = [[0.13, -0.21] .+ G .* t for t in range(0, 1, 9)]
                expected = sort(angle.(cis.([2π * dot(G, r) for r in rs])))
                @test all(phase_isapprox.(wilsonloop(ptbm′, ks, 1:N), expected))
                N == 1 && @test phase_isapprox(wilsonloop(ptbm′, ks, 1), only(expected))
            end
        end
    end

    @testset "Kramers degeneracy on time-reversal invariant lines (Kane-Mele model)" begin
        brs′ = bandreps(17, Val(2); spinful = Val(true)) # p6mm, as in docs/src/spinful.md
        tbm′ = tb_hamiltonian(@composite(brs′[3]), [[0, 0], [1, 0]])[1:6]
        k2s = range(-0.5, 0.5, 41)
        # Kane-Mele, without & with Rashba SOC; lowest Kramers pair (bands 1:2) is isolated
        for ptbm′ in (tbm′([0, 0, 1, 0, 0, 0.05]), tbm′([0, 0.05, 1, 0, 0, 0.05]))
            for k1 in (0.0, 0.5) # time-reversal invariant lines: Kramers-degenerate phases
                θs = wilsonloop(ptbm′, [[k1, k2] for k2 in k2s], 1:2)
                @test phase_isapprox(θs[1], θs[2])
            end
            θs = wilsonloop(ptbm′, [[0.25, k2] for k2 in k2s], 1:2) # generic line: split
            @test !phase_isapprox(θs[1], θs[2]; atol = 1e-2)
        end
    end

    # an isolated band with a Wannier function centered at q gives a Wilson loop phase
    # 2πG⋅q along a loop of winding G; since Convention 1 includes orbital positions in the
    # Bloch basis, q is absolute (not relative to the orbitals). We identify q independently
    # via the band's symmetry content: as the Wyckoff position of the matching EBR
    @testset "Wannier-center interpretation" begin
        function test_wannier_center(ptbm, brs, Gs, k0s, ebr_label)
            br = lowest_band_ebr(ptbm, brs)
            @test string(br) == ebr_label
            q = constant(position(br))
            for G in Gs, k0 in k0s
                ks = [k0 .+ G .* t for t in range(0, 1, 41)]
                @test phase_isapprox(wilsonloop(ptbm, ks, 1), 2π * dot(G, q))
            end
        end

        @testset "1D: Zak phase" begin
            brs1 = bandreps(2, Val(1)) # p1̄; 1a = [0], 1b = [1/2]
            Gs, k0s = ([1], [2], [-1]), ([0.0], [0.37])
            # s orbitals at 1b & 1a (terms: ε_1b, t_1b, ε_1a, t_1a, t_ab): band at the lower site
            tbm_ab = tb_hamiltonian(@composite(brs1[1] + brs1[3]), [[0], [1]])
            test_wannier_center(tbm_ab([-1, 0.2, 1, 0.1, 0.3]), brs1, Gs, k0s, "(1b|A′)")
            test_wannier_center(tbm_ab([1, 0.2, -1, 0.1, 0.3]), brs1, Gs, k0s, "(1a|A′)")
            # s & p orbitals at 1b (terms: ε_s, t_ss, ε_p, t_pp, t_sp): band at 1b if trivial,
            # but at 1a if inverted, i.e., an obstructed atomic limit (no orbitals at 1a)
            tbm_sp = tb_hamiltonian(@composite(brs1[1] + brs1[2]), [[0], [1]])
            test_wannier_center(tbm_sp([-1, 0.2, 1, -0.2, 0.3]), brs1, Gs, k0s, "(1b|A′)")
            test_wannier_center(tbm_sp([0, -0.5, 0, 0.5, 0.3]), brs1, Gs, k0s, "(1a|A′)")
        end

        # in 2D, the phase of a loop along G at transverse offset k⊥ is a hybrid Wannier
        # center, which generally varies with k⊥; it is pinned to G⋅q on loops through
        # C₂-invariant momenta (`trims`) and, if the band is isolated over the entire BZ
        # and C₂T is a symmetry (as here, with time-reversal), on every loop (`generic`)
        @testset "2D: hybrid Wannier centers (p2)" begin
            brs2 = bandreps(2, Val(2)) # 1a = [0,0], 1b = [0,1/2], 1c = [1/2,0], 1d = [1/2,1/2]
            Gs = ([1,0], [0,1], [1,1], [2,1])
            trims = ([0.0,0], [0.5,0], [0,0.5], [0.5,0.5])
            generic = ([0.0,0.31], [0.23,0.0], [0.17,-0.41])
            # A orbitals at all of 1d, 1c, 1b, 1a, in an atomic limit: band at the lowest site
            tbm4 = tb_hamiltonian(CompositeBandRep([1,0,1,0,1,0,1,0], brs2), [[0,0], [1,0], [0,1]])
            onsite_idxs = [1, 4, 7, 10] # on-site terms of 1d, 1c, 1b, 1a
            for (j, label) in enumerate(("(1d|A)", "(1c|A)", "(1b|A)", "(1a|A)"))
                cs = [0.1*cospi(0.73*k) for k in 1:length(tbm4)]
                cs[onsite_idxs] .= [0.5, 2.0, 3.5, 5.0]
                cs[onsite_idxs[j]] = -1.0 # make site `j` the lowest-energy site
                test_wannier_center(tbm4(cs), brs2, Gs, (trims..., generic...), label)
            end
            # A & B orbitals at 1a, inverted along δ = [1,0] (terms 2 & 6): band at 1c, i.e.,
            # an obstructed atomic limit (small asymmetric terms avoid C₂T-protected Dirac
            # points, which a more symmetric choice of hoppings produces)
            tbm2 = tb_hamiltonian(CompositeBandRep([0,0,0,0,0,0,1,1], brs2),
                                  [[0,0], [1,0], [0,1], [1,1]])
            ptbm2 = tbm2([0.0, -0.45, -0.03, 0.02, 0.0, 0.46, 0.03, 0.01, 0.3, 0.25, 0.1])
            @test lowest_gap(ptbm2) > 0.05
            test_wannier_center(ptbm2, brs2, Gs, (trims..., generic...), "(1c|A)")
        end
    end

    @testset "Input forms" begin
        ptbm = haldane_model(1.0, 0.5, 0.3, π/2)
        k2s = range(-0.5, 0.5, 21)
        ks = [[0.1, k2] for k2 in k2s]
        γ = wilsonloop(ptbm, ks, 1)
        θs = wilsonloop(ptbm, ks, 1:2)
        # k-points as vectors, tuples, `ReciprocalPoint`s, or a lazy iterator
        @test wilsonloop(ptbm, [(0.1, k2) for k2 in k2s], 1) ≈ γ
        @test wilsonloop(ptbm, [ReciprocalPoint(0.1, k2) for k2 in k2s], 1) ≈ γ
        @test wilsonloop(ptbm, ([0.1, k2] for k2 in k2s), 1) ≈ γ
        # integer vs. single-element band selections; order of band selection is irrelevant
        @test wilsonloop(ptbm, ks, [1]) == [γ]
        @test wilsonloop(ptbm, ks, [2, 1]) ≈ wilsonloop(ptbm, ks, [1, 2]) ≈ θs
        # the starting point of the loop is irrelevant (same points, cyclically shifted)
        ks_shifted = [[0.1, k2] for k2 in range(-0.3, 0.7, 21)]
        @test phase_isapprox(wilsonloop(ptbm, ks_shifted, 1), γ)
        @test all(phase_isapprox.(wilsonloop(ptbm, ks_shifted, 1:2), θs))
    end

    @testset "Invalid input" begin
        ptbm = haldane_model(1.0, 0.5, 0.3, π/2)
        ks = [[0.1, k2] for k2 in range(-0.5, 0.5, 21)]
        @test_throws ErrorException wilsonloop(ptbm, [[0.1, k2] for k2 in range(0, 0.5, 21)], 1) # open loop
        @test_throws ErrorException wilsonloop(ptbm, [[0.1, k2, 0] for k2 in range(-0.5, 0.5, 21)], 1) # wrong dim.
        @test_throws ErrorException wilsonloop(ptbm, [[0.1, 0.0]], 1) # too few k-points
        @test_throws ErrorException wilsonloop(ptbm, ks, 3)      # band out of range
        @test_throws ErrorException wilsonloop(ptbm, ks, Int[])  # no bands
        @test_throws ErrorException wilsonloop(ptbm, ks, [1, 1]) # repeated bands
        tbm_nh = tb_hamiltonian(cbr, [[0,0]], Val(NONHERMITIAN))
        ptbm_nh = tbm_nh([0.3*cospi(0.73*k) for k in 1:length(tbm_nh)])
        @test_throws ErrorException wilsonloop(ptbm_nh, ks, 1)   # non-Hermitian
    end
end
