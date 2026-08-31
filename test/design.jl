using Test
using SymmetricTightBinding
using SymmetricTightBinding: solve, FAILED_ISOLATION_LOSS, _reduced_distance
using Crystalline
using Optim
using Random

# Tests for the objective-driven design tools (`src/design.jl`) and the `isolate_irrep`
# driver in the Optim.jl extension. The workhorse model is a deliberately tiny 2D one — a
# p6mm (1a|E₁) ⊕ (1a|A₁) pair, i.e., 3 bands and 6 terms, whose Γ-point content is Γ₁+Γ₅ —
# so that the target multiplet (the 2-dimensional Γ₅) must split one band below and one
# above its energy, with the Γ₁ singlet joining one side.

@testset "Design (irrep isolation)" begin
    brs = calc_bandreps(17, Val(2))
    cbr = @composite brs[12] + brs[8]                  # (1a|E₁) ⊕ (1a|A₁): Γ₁+Γ₅, 3 bands
    tbm = tb_hamiltonian(cbr, [[0, 0], [1, 0]])        # on-site + nearest neighbour
    Nᶜ = length(tbm)
    csv = [[0.3cospi(0.73(k + 7j)) for k in 1:Nᶜ] for j in 0:4] # deterministic amplitudes
    # a few generic (non-high-symmetry) k-points, for the gradient checks: at a symmetry-
    # enforced degeneracy the sorted band energies are only one-sidedly differentiable, so
    # central differences would not converge to the analytic gradient there
    ks_gen = [[0.21, 0.13], [0.37, -0.29], [-0.11, 0.42], [0.44, 0.05]]

    # -------------------------------------------------------------------------------------
    @testset "uniform_kmesh" begin
        ks = uniform_kmesh(Val(2), 4)
        @test length(ks) == 16
        @test count(k -> iszero(k), ks) == 1                     # Γ is sampled, exactly once
        @test all(k -> all(x -> -0.5 ≤ x < 0.5, k), ks)          # inside the reduced BZ
        @test length(uniform_kmesh(Val(3), 3)) == 27
        @test_throws ErrorException uniform_kmesh(Val(2), 0)
    end

    @testset "_reduced_distance (minimal image)" begin
        Γ = [0.0, 0.0]
        @test _reduced_distance([0.3, 0.0], Γ) ≈ 0.3
        @test _reduced_distance([-0.3, 0.0], Γ) ≈ 0.3
        @test _reduced_distance([0.9, 0.0], Γ) ≈ 0.1              # wraps to -0.1
        @test _reduced_distance([0.5, 0.5], Γ) ≈ sqrt(0.5)
        @test _reduced_distance([0.2, 0.1], [0.2, 0.1]) ≈ 0.0
    end

    # -------------------------------------------------------------------------------------
    @testset "IrrepTarget" begin
        target = IrrepTarget(tbm, "Γ₅")
        @test target.d == 2
        @test iszero(target.k)
        @test length(target.Ms) == length(target.χ) == 12         # |6mm| = 12 operations
        @test_throws ErrorException IrrepTarget(tbm, "Γ₉")        # no such irrep
        @test_throws ErrorException IrrepTarget(tbm, "Γ₆")        # absent from `cbr`'s content
    end

    # -------------------------------------------------------------------------------------
    @testset "locate_multiplet agrees with `collect_irrep_annotations`" begin
        target = IrrepTarget(tbm, "Γ₅")
        for cs in csv
            ptbm = tbm(cs)
            Es, vs = solve(ptbm, target.k; bloch_phase = Val(false))
            annotations = collect_irrep_annotations(ptbm)["Γ"]
            expected = first(only(filter(p -> last(p) == "Γ₅", annotations)))
            @test locate_multiplet(target, Es, vs) == expected
        end
    end

    # -------------------------------------------------------------------------------------
    @testset "IrrepIsolationObjective" begin
        obj = IrrepIsolationObjective(tbm, "Γ₅"; ks = ks_gen)

        @testset "gradient matches finite differences" begin
            # (also with an ℓ₁ penalty active: `csv`'s amplitudes are all nonzero, so |c| is
            # differentiable at each of them)
            for obj′ in (obj, IrrepIsolationObjective(tbm, "Γ₅"; ks = ks_gen, lasso = 0.3))
                for cs in csv
                    G = zeros(Nᶜ)
                    obj′(1.0, G, nothing, cs)
                    Gᶠᵈ = map(1:Nᶜ) do i
                        δ = zeros(Nᶜ); δ[i] = 1e-6
                        (obj′(1.0, nothing, nothing, cs + δ) -
                         obj′(1.0, nothing, nothing, cs - δ)) / 2e-6
                    end
                    @test G ≈ Gᶠᵈ rtol=1e-5
                end
            end
        end

        @testset "lasso penalty" begin
            objₗ = IrrepIsolationObjective(tbm, "Γ₅"; ks = ks_gen, lasso = 0.3)
            for cs in csv
                # the penalized loss exceeds the unpenalized one by exactly λ₁‖c‖₁/(Nᶜσ),
                # with σ the spectrum's RMS spread (available from the report)
                σ = isolation_report(tbm(cs), "Γ₅"; ks = ks_gen).σ
                Δ = objₗ(1.0, nothing, nothing, cs) - obj(1.0, nothing, nothing, cs)
                @test Δ ≈ 0.3 * sum(abs, cs) / (Nᶜ * σ)
            end
            # and the penalty is scale-invariant, like the rest of the loss
            objₗ′ = IrrepIsolationObjective(tbm, "Γ₅"; ks = ks_gen, lasso = 0.3,
                                            λ_scale = 0, λ_shift = 0)
            @test objₗ′(1.0, nothing, nothing, 2.5csv[1]) ≈
                  objₗ′(1.0, nothing, nothing, csv[1])
        end

        @testset "scale & offset invariance of the (unpenalized) loss" begin
            # with the gauge penalties off, the loss must be invariant under an overall
            # rescaling of the amplitudes and under an overall energy shift (here: terms 1 &
            # 3, the on-site energies of the two EBRs, which sum to the identity)
            obj′ = IrrepIsolationObjective(tbm, "Γ₅"; ks = ks_gen, λ_scale = 0, λ_shift = 0)
            cs = csv[1]
            L = obj′(1.0, nothing, nothing, cs)
            @test obj′(1.0, nothing, nothing, 2.5cs) ≈ L
            shift = zeros(Nᶜ); shift[1] = 1.0; shift[3] = 1.0 # ↔ H(k) → H(k) + 1
            @test obj′(1.0, nothing, nothing, cs + shift) ≈ L
        end

        @testset "diagnostics & failure sentinel" begin
            obj(1.0, nothing, nothing, csv[1])
            @test length(obj.bands) == 2
            @test obj.split == 1        # the only proper split of a 2-dimensional irrep
            @test isfinite(obj.relgap)

            # a vanishing Hamiltonian leaves every band degenerate, so the Γ₅ multiplet
            # cannot be told apart from the Γ₁ band: the objective must bail out
            @test obj(1.0, nothing, nothing, zeros(Nᶜ)) == FAILED_ISOLATION_LOSS
            @test isempty(obj.bands)

            @test_throws ErrorException obj(1.0, nothing, zeros(Nᶜ, Nᶜ), csv[1])
        end

        @testset "constraints inside the excluded neighbourhood" begin
            # `kexclude` exempts near-`k₀` points from the *multiplet-splitting* constraints
            # only: the bands outside the multiplet stay a finite distance from `Et` there,
            # and must still be pushed away — so such points are constrained, just less
            near, far = [0.03, 0.0], [0.31, 0.17]
            obj′ = IrrepIsolationObjective(tbm, "Γ₅"; ks = [near, far], kexclude = 0.1)
            obj′(1.0, nothing, nothing, csv[1])         # populates `bands`/`split`
            nc = SymmetricTightBinding._constraints!(obj′, obj′.bands, obj′.split)
            per_κ = [count(c -> first(c) == κ, view(obj′.cs_buf, 1:nc)) for κ in 1:3]
            @test per_κ[1] == per_κ[2]                  # k₀ & `near`: outside bands only
            @test per_κ[3] == per_κ[1] + 2              # `far`: + the two split constraints
        end

        @testset "input validation" begin
            # a model whose every band belongs to the multiplet poses no isolation problem
            @test_throws ErrorException IrrepIsolationObjective(tb_hamiltonian(
                (@composite brs[12]), [[0, 0], [1, 0]]), "Γ₅")
            # nor does a one-dimensional irrep, which cannot be split across its own energy
            @test_throws ErrorException IrrepIsolationObjective(tbm, "Γ₁")
            @test_throws ErrorException IrrepIsolationObjective(tbm, "Γ₅"; kexclude = 0)
            @test_throws ErrorException IrrepIsolationObjective(tbm, "Γ₅"; split = 2)
            # and neither does a k-sampling that the exclusion zone empties out entirely
            @test_throws ErrorException IrrepIsolationObjective(tbm, "Γ₅";
                                            ks = [[0.02, 0.0]], kexclude = 0.1)
        end
    end

    # -------------------------------------------------------------------------------------
    @testset "isolate_irrep & isolation_report" begin
        Random.seed!(1)
        # NB: `Nk` divisible by 6, so that the mesh samples K = (⅓, ⅓) — where the K₃ irrep
        # forces a two-fold degeneracy that rules out one of the two fillings outright. A
        # mesh missing K (e.g. `Nk = 8`) yields models that only *look* isolated.
        ptbm = isolate_irrep(tbm, "Γ₅"; Nk = 12, max_multistarts = 60)

        # the outcome is verified independently, on a mesh 2.5× finer than the one optimized
        # over, and against the irrep labels that `collect_irrep_annotations` assigns
        report = isolation_report(ptbm, "Γ₅"; Nk = 30)
        @test report.isolated
        @test report.gap > 0 && report.relgap > 0
        @test report.split == 1
        @test report.bands == first(only(filter(p -> last(p) == "Γ₅",
                                                collect_irrep_annotations(ptbm)["Γ"])))

        # the isolation condition itself: away from Γ, the multiplet's own bands must sit on
        # either side of `Et`, one below (band ν) and one above (band ν+1) — and hence so
        # must every other band, by sorting
        ν = first(report.bands) + report.split - 1
        for k in uniform_kmesh(Val(2), 30)
            _reduced_distance(k, [0.0, 0.0]) ≥ 0.1 || continue
            Es = spectrum_single_k(ptbm, k)
            @test Es[ν] < report.Et < Es[ν+1]
        end

        # a model whose multiplet cannot be located at all is reported as not isolated
        report′ = isolation_report(tbm(zeros(Nᶜ)), "Γ₅"; Nk = 5)
        @test !report′.isolated && isempty(report′.bands)
    end
end
