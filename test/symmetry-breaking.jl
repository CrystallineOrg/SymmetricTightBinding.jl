using Test
using SymmetricTightBinding
using SymmetricTightBinding: _subduced_complement, _issubgroup, _subduction_groups
using Crystalline

@testset "Symmetry breaking" begin
    @testset "2D example from docs" begin
        D = 2
        brs = bandreps(11, Val(D); timereversal = true)
        cbr = @composite brs[1] # (2c|A₁)
        Rs = [[0,0], [1,0]]
        tbm = tb_hamiltonian(cbr, Rs)
        
        Δtbm_C4  = subduced_complement(tbm, Rs, 6)                        # break C₄
        Δtbm_m   = subduced_complement(tbm, Rs, 10)                       # break mirror
        Δtbm_tr  = subduced_complement(tbm, Rs, 11; timereversal = false) # break TR
        Δtbm_mtr = subduced_complement(tbm, Rs, 10; timereversal = false) # break both
        @test length(Δtbm_C4) == 4
        @test length(Δtbm_m) == 1
        @test length(Δtbm_tr) == 1
        @test length(Δtbm_mtr) == 4
        
        # adding more orbits, we find another mirror-breaking term, but no TR-breaking one
        Rs_big = [[0,0], [1,0], [1,1]]
        tbm_big = tb_hamiltonian(cbr, Rs_big)
        Δtbm_big_m   = subduced_complement(tbm_big, Rs_big, 10)                       # break mirror
        Δtbm_big_tr  = subduced_complement(tbm_big, Rs_big, 11; timereversal = false) # break TR
        Δtbm_big_mtr = subduced_complement(tbm_big, Rs_big, 10; timereversal = false) # break both
        @test length(Δtbm_big_m) == 2
        @test length(Δtbm_big_tr) == 1
        @test length(Δtbm_big_mtr) == 5
        @test issubset(Δtbm_big_m, Δtbm_big_mtr) # must subset eachother

        # breaking mirror and TR symmetry together should give the same basis as starting
        # directly with plane group p4 (#10) and breaking TR from the get-go
        brs10 = bandreps(10, Val(D); timereversal=false)
        cbr10 = @composite brs10[1] # (2c|A) (unlike #11: two distinct M irreps, M₃ & M₄)

        tbm10 = tb_hamiltonian(cbr10, [[0,0], [1,0], [0,1]]) # (⋆)
        @test length(tbm10) == length(tbm) + length(Δtbm_mtr)
        # (⋆): must add [0,1] also, cf. diagonally-directed hopping term involving both 
        # [1,0] and [0,1] in C₄ setting

        tbm10_big = tb_hamiltonian(cbr10, Rs_big)
        @test length(tbm10_big) == length(tbm_big) + length(Δtbm_big_mtr)
    end

    @testset "3D example" begin
        # this is not a well-thought out example, but just added to test that things work
        # without erroring
        brs = bandreps(16, Val(3); timereversal = true)
        cbr = @composite brs[1] + brs[end] #  (1h|A) + (1a|B₂) (2 bands)
        Rs = [[0,0,0],]
        tbm = tb_hamiltonian(cbr, Rs)

        @test length(subduced_complement(tbm, Rs, 3)) == 1
    end

    @testset "3D example, body-centered (I) lattice" begin
        brs = bandreps(121, Val(3); timereversal = true)
        cbr = @composite brs[1] + brs[end-1] # (4d|A) + (2a|A₂) (3 bands)
        Rs = [[0,0,0],]
        tbm = tb_hamiltonian(cbr, Rs)

        Δtbm = subduced_complement(tbm, Rs, 82)
        @test length(Δtbm) == 2

        # completeness: `(4d|A)` of ⋕121 splits into `2c ⊕ 2d` in ⋕82, and `(2a|A₂)` maps
        # to one of `(2a|A/B)`; every such 3-band composite of ⋕82 has 6 = 4 + 2 terms
        brs82 = bandreps(82, Val(3); timereversal = true)
        for (i, j, k) in Iterators.product((4, 5), (1, 2), (10, 11)) # 2c, 2d, 2a
            cbr82 = CompositeBandRep([n ∈ (i,j,k) ? 1 : 0 for n in eachindex(brs82)], brs82)
            @test length(tb_hamiltonian(cbr82, Rs)) == length(tbm) + length(Δtbm)
        end

        # having added the complement, there is nothing further to find in ⋕82
        @test length(subduced_complement(vcat(tbm, Δtbm), Rs, 82)) == 0
    end

    @testset "composite band representations at a shared Wyckoff position" begin
        # terms belonging to *different* blocks can carry equal (`==`) hopping orbits, if
        # the associated band representations sit at the same Wyckoff position; such terms
        # must not be grouped together, since they do not share a coefficient basis
        brs = bandreps(47, Val(3); timereversal = true) # P4/mmm
        cbr = @composite brs[57] + brs[60] # (1a|Ag) + (1a|B₁ᵤ)

        Rs₀ = [[0,0,0]]
        tbm = tb_hamiltonian(cbr, Rs₀) # on-site only: blocks (1,1) & (2,2), both δ=0
        @test length(tbm) == 2
        # subducing to G itself, with unchanged time-reversal, must give nothing new
        @test length(subduced_complement(tbm, Rs₀, 47)) == 0
        @test length(subduced_complement(tbm, Rs₀, 47; timereversal = false)) == 0

        Rs = [[0,0,0], [1,0,0], [0,1,0], [0,0,1]]
        tbm_nn = tb_hamiltonian(cbr, Rs)
        @test length(tbm_nn) == 9
        @test length(subduced_complement(tbm_nn, Rs, 47)) == 0
        Δtbm_nn = subduced_complement(tbm_nn, Rs, 47; timereversal = false)
        @test length(Δtbm_nn) == 1

        # the complement must be *complete*: breaking time-reversal in G should give the
        # same number of terms as building the model without time-reversal from the start
        brs′ = bandreps(47, Val(3); timereversal = false)
        cbr′ = @composite brs′[57] + brs′[60]
        # ⋕57 & ⋕60 index the same band representations with and without time-reversal, but
        # that is a property of `bandreps`' ordering rather than something we control
        @test string.((brs[57], brs[60])) == ("(1a|Ag)", "(1a|B₁ᵤ)")
        @test string.((brs′[57], brs′[60])) == ("(1a|Ag)", "(1a|B₁ᵤ)")
        @test length(tb_hamiltonian(cbr′, Rs)) == length(tbm_nn) + length(Δtbm_nn)

        # the extended model must still be Hermitian
        tbm′ = vcat(tbm_nn, Δtbm_nn)
        ptbm′ = tbm′(rand(length(tbm′)))
        for k in (ReciprocalPoint(0.1, 0.2, 0.3), ReciprocalPoint(0.5, 0.0, 0.25))
            @test ptbm′(k) ≈ ptbm′(k)' # NB: `ptbm(k)` returns a reused buffer
        end

        # repeated band representations: blocks (1,1), (2,2), and (1,2) all share orbits
        cbr2 = @composite brs[57] + brs[57] # 2 × (1a|Ag)
        Rs₁ = [[0,0,0], [1,0,0]]
        tbm2 = tb_hamiltonian(cbr2, Rs₁)
        @test length(tbm2) == 6
        @test length(subduced_complement(tbm2, Rs₁, 47)) == 0
        @test length(subduced_complement(tbm2, Rs₁, 47; timereversal = false)) == 2
    end

    @testset "term grouping" begin
        # terms sharing a coefficient basis must be grouped together regardless of whether
        # they appear contiguously in the model (models may be built by `vcat` or indexing)
        brs2d = bandreps(11, Val(2); timereversal = true)
        Rs = [[0,0], [1,0]]
        tbm = tb_hamiltonian((@composite brs2d[1]), Rs)
        groupidxs(tbm, Rs) = [g.idxs for g in _subduction_groups(tbm, Rs)]
        @test groupidxs(tbm, Rs) == [[1], [2], [3, 4], [5]]

        p = [3, 1, 2, 4, 5] # splits the {3,4} group apart
        @test groupidxs(tbm[p], Rs) == [[2], [3], [1, 4], [5]]
        for (sgnumᴴ, timereversal) in ((10, true), (6, true), (11, false), (10, false))
            @test length(subduced_complement(tbm[p], Rs, sgnumᴴ; timereversal)) ==
                  length(subduced_complement(tbm, Rs, sgnumᴴ; timereversal))
        end

        # equal-orbit terms from distinct blocks must stay in distinct groups, also when
        # they are adjacent
        brs = bandreps(47, Val(3); timereversal = true)
        Rs₁ = [[0,0,0], [1,0,0]]
        tbm2 = tb_hamiltonian((@composite brs[57] + brs[57]), Rs₁)
        q = [1, 3, 5, 2, 4, 6] # interleave, so that equal orbits become adjacent
        @test [t.block_ij for t in tbm2.terms] ==
              [(1,1), (1,1), (2,2), (2,2), (1,2), (1,2)]
        @test groupidxs(tbm2[q], Rs₁) == [[1], [4], [2], [5], [3], [6]]
        @test length(subduced_complement(tbm2[q], Rs₁, 47; timereversal = false)) ==
              length(subduced_complement(tbm2,    Rs₁, 47; timereversal = false))
    end

    @testset "centered lattices" begin
        # the subgroup generators must be converted to the primitive setting before the
        # constraints are imposed; if not, centered lattices error out in
        # `site_induced_sgrep` (which compares against primitivized site groups)
        brs = bandreps(12, Val(3); timereversal = true) # C2/m (C-centered)
        cbr = @composite brs[1] # (4f|Ag)
        Rs = [[0,0,0], [1,0,0]]
        tbm = tb_hamiltonian(cbr, Rs)
        @test length(tbm) == 6

        # subducing to G itself, with unchanged time-reversal, must give nothing new
        @test length(subduced_complement(tbm, Rs, 12)) == 0

        @test length(subduced_complement(tbm, Rs, 12; timereversal = false)) == 1 # break TR
        @test length(subduced_complement(tbm, Rs, 5)) == 2                       # break mirror
        @test length(subduced_complement(tbm, Rs, 8)) == 2                       # break C₂ & -1
    end

    @testset "orbits that are empty in the parent group (issue #117)" begin
        # a (block, orbit) pair on which *every* coefficient is forbidden in G carries no
        # term, and so cannot be found without `Rs` - even though the symmetry reduction
        # may be exactly what allows it
        Rs = [[0,0,0], [1,0,0], [0,1,0], [0,0,1]]
        brs = bandreps(47, Val(3); timereversal = true) # P4/mmm
        cbr = @composite brs[57] + brs[60] # (1a|Ag) + (1a|B₁ᵤ): s & p_z on a shared site
        tbm = tb_hamiltonian(cbr, Rs)
        @test length(tbm) == 9

        # the (1,2) block is forbidden on the x and y bonds in P4/mmm; dropping m_z (while
        # keeping inversion and TR) frees the bond along the retained 2-fold axis
        gensᴴ = [S"x,-y,-z", S"-x,-y,-z"] # 2ₓ & -1, i.e. 2/m with unique axis a
        Δtbm = _subduced_complement(tbm, Rs, gensᴴ)
        @test length(Δtbm) == 1
        @test only(Δtbm).block_ij == (1, 2)
        @test length(subduced_complement(tbm, Rs, 10)) == 1      # ⋕10 = P2/m

        # completeness: the direct P2/m model over the same range has exactly one more term
        brs10 = bandreps(10, Val(3); timereversal = true)
        @test string.((brs10[29], brs10[32])) == ("(1a|Ag)", "(1a|Bᵤ)")
        tbm10 = tb_hamiltonian((@composite brs10[29] + brs10[32]), Rs)
        @test length(tbm10) == length(tbm) + length(Δtbm)

        # the extended model must still be Hermitian, and have nothing further to give
        tbm′ = vcat(tbm, Δtbm)
        ptbm′ = tbm′(rand(length(tbm′)))
        for k in (ReciprocalPoint(0.1, 0.2, 0.3), ReciprocalPoint(0.5, 0.0, 0.25))
            @test ptbm′(k) ≈ ptbm′(k)' # NB: `ptbm(k)` returns a reused buffer
        end
        @test length(subduced_complement(tbm′, Rs, 10)) == 0

        # the pairs over `Rs` include some that carry no term of `tbm`, and every term of
        # `tbm` belongs to exactly one pair
        groups = _subduction_groups(tbm, Rs)
        @test count(g -> isempty(g.idxs), groups) > 0
        @test sort(reduce(vcat, g.idxs for g in groups)) == 1:length(tbm)
    end

    @testset "hopping range `Rs`" begin
        brs = bandreps(11, Val(2); timereversal = true) # p4mm
        cbr = @composite brs[1] # (2c|A₁)
        Rs = [[0,0], [1,0]]
        tbm = tb_hamiltonian(cbr, Rs)

        @test length(subduced_complement(tbm, Rs, 11)) == 0 # subducing to G itself

        # `Rs` alone sets the (block, orbit) pairs that are searched: a narrower `Rs` gives
        # exactly the complement of the correspondingly narrower model, and an empty one
        # gives nothing
        Rs′ = [[0,0]]
        tbm_narrow = tb_hamiltonian(cbr, Rs′)
        for (sgnumᴴ, timereversal) in ((11, true), (10, true), (10, false), (6, true))
            @test subduced_complement(tbm, Rs′, sgnumᴴ; timereversal).terms ==
                  subduced_complement(tbm_narrow, Rs′, sgnumᴴ; timereversal).terms
        end
        @test isempty(subduced_complement(tbm, Vector{Int}[], 10))

        # the tests below pin what happens for unsupported usage - a sub-selected `tbm`, or
        # a wider `Rs` than the model's own - so that the behavior is at least deterministic

        # a sub-selected `tbm` gives the complement relative to the sub-selection, so terms
        # dropped by it return as "new" even when subducing to G itself
        @test length(subduced_complement(tbm[[1,2,3,5]], Rs, 11)) == 1 # dropped term 4
        @test length(subduced_complement(tbm[1:4], Rs, 11)) == 1       # dropped term 5

        # a wider `Rs` returns the longer-range terms too, whether or not G allows them
        Rs_big = [[0,0], [1,0], [1,1]]
        @test length(subduced_complement(tbm, Rs_big, 11)) ==
              length(tb_hamiltonian(cbr, Rs_big)) - length(tbm)

        # completeness against a directly-built subgroup model: the 2c orbit of p4mm splits
        # into 1c ⊕ 1b of p2mm. The orbits are those of p4mm, and reach further than `Rs`,
        # so the ⋕6 model must be built over a range covering the same hopping vectors
        brs6 = bandreps(6, Val(2); timereversal = true)
        cbr6 = CompositeBandRep([n ∈ (5, 9) ? 1 : 0 for n in eachindex(brs6)], brs6)
        @test string(cbr6) == "(1c|A₁) + (1b|A₁)"
        @test length(tb_hamiltonian(cbr6, [[0,0], [1,0], [0,1], [-1,0]])) ==
              length(tbm) + length(subduced_complement(tbm, Rs, 6))

        # non-Hermitian models iterate over all blocks, not just the upper-triangular ones
        tbm_nh = tb_hamiltonian(cbr, Rs, Val(NONHERMITIAN))
        @test length(subduced_complement(tbm_nh, Rs, 11)) == 0
        @test length(subduced_complement(tbm_nh, Rs, 10)) == 3
        @test length(subduced_complement(tbm_nh, Rs, 6)) == 6
    end

    @testset "subgroup precondition" begin
        # `_subduced_complement` asserts that `gensᴴ` generate a subgroup of G, given in G's
        # conventional setting; unreachable via `subduced_complement`, but check it can fire
        brs2d = bandreps(11, Val(2); timereversal = true) # p4mm
        Rs = [[0,0], [1,0]]
        tbm = tb_hamiltonian((@composite brs2d[1]), Rs)

        @test _issubgroup([S"y,x"], 11)     # mₓᵧ ∈ p4mm
        @test !_issubgroup([S"-y,x+y"], 11) # 6⁺ ∉ p4mm
        @test length(_subduced_complement(tbm, Rs, [S"y,x"])) isa Int
        @test_throws AssertionError _subduced_complement(tbm, Rs, [S"-y,x+y"])

        # for double group operations, the SU(2) elements are checked as well
        g = S"x,-y,-z" # 2₁₀₀ ∈ Pmmm (⋕47)
        @test _issubgroup([DSymOperation{3}(g, SU2(g, 47))], 47)
        @test !_issubgroup([DSymOperation{3}(g, one(SU2))], 47)
    end

    @testset "spinful models" begin
        # completeness: breaking time-reversal in P-1 gives as many terms as building the
        # model without time-reversal from the start; the Kramers pair (1h|AᵤˢAᵤˢ) splits
        # into two copies of (1h|Aᵤˢ) without time-reversal
        Rs1 = [[0,0,0], [1,0,0], [0,1,0], [0,0,1]]
        brs = bandreps(2; spinful = Val(true), timereversal = true)
        tbm = tb_hamiltonian((@composite brs[1]), Rs1) # (1h|AᵤˢAᵤˢ)
        Δtbm = subduced_complement(tbm, Rs1, 2; timereversal = false)
        brs′ = bandreps(2; spinful = Val(true), timereversal = false)
        @test string(brs[1]) == "(1h|AᵤˢAᵤˢ)" && string(brs′[1]) == "(1h|Aᵤˢ)"
        @test length(tb_hamiltonian((@composite 2brs′[1]), Rs1)) ==
              length(tbm) + length(Δtbm)

        # spinful models require double group operations (and vice versa)
        @test_throws "requires a double group operation" _subduced_complement(tbm, Rs1, [S"-x,-y,-z"])

        # the full model `vcat(tbm, Δtbm)` is symmetric under exactly the operations of H's
        # double group (among those of G's), and remains time-reversal symmetric
        ks = [ReciprocalPoint(0.13, 0.27, 0.19), ReciprocalPoint(0.31, -0.07, 0.44)]
        Rs0 = [[0,0,0], [1,0,0]]
        for (sgnumᴳ, sgnumᴴ, idxs, Rs) in (
                (191, 183, [1],    Rs0),        # P6/mmm → P6mm: hexagonal
                (225, 216, [1],    [[0,0,0]]),  # Fm-3m → F-43m: F centring (24d: on-site
                                                #   only, since longer range is slow)
                (47,  25,  [1, 2], Rs0))        # Pmmm → Pmm2: composite
            brs = bandreps(sgnumᴳ; spinful = Val(true), timereversal = true)
            cbr = CompositeBandRep([count(==(n), idxs) for n in eachindex(brs)], brs)
            tbm = tb_hamiltonian(cbr, Rs)
            Δtbm = subduced_complement(tbm, Rs, sgnumᴴ)
            @test length(Δtbm) > 0
            tbm′ = vcat(tbm, Δtbm)
            ptbm′ = tbm′([0.3*cospi(0.73*n) for n in 1:length(tbm′)])
            opsᴳ = primitivize(spacegroup(sgnumᴳ, Val(3); spinful = Val(true)))
            nopsᴴ = length(primitivize(spacegroup(sgnumᴴ, Val(3); spinful = Val(true))))
            Γ = SymmetricTightBinding.site_induced_timereversal_unitary(tbm′.cbr)
            preserved = map(opsᴳ) do g
                all(ks) do k
                    Hk = copy(ptbm′(k))
                    D = site_induced_sgrep(ptbm′, g)(k)
                    isapprox(ptbm′(g * k), D * Hk * D'; atol = 1e-10)
                end
            end
            @test count(preserved) == nopsᴴ
            @test all(ks) do k
                Hk = copy(ptbm′(k))
                isapprox(ptbm′(-k), Γ * conj(Hk) * Γ'; atol = 1e-10)
            end
            # having added the complement, there is nothing further to find in H
            @test length(subduced_complement(tbm′, Rs, sgnumᴴ)) == 0
        end
    end
end
