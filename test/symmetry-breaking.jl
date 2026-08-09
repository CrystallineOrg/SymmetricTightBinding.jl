using Test
using SymmetricTightBinding
using SymmetricTightBinding: _group_terms_by_block_and_orbit, _subduced_complement,
                             _issubgroup, _subduction_groups
using Crystalline

@testset "Symmetry breaking" begin
    @testset "2D example from docs" begin
        D = 2
        brs = calc_bandreps(11, Val(D); timereversal = true)
        cbr = @composite brs[1] # (2c|A₁)
        tbm = tb_hamiltonian(cbr, [[0,0], [1,0]])
        
        Δtbm_m   = subduced_complement(tbm, 10)                       # break mirror
        Δtbm_tr  = subduced_complement(tbm, 11; timereversal = false) # break TR
        Δtbm_mtr = subduced_complement(tbm, 10; timereversal = false) # break both
        @test length(Δtbm_m) == 1
        @test length(Δtbm_tr) == 1
        @test length(Δtbm_mtr) == 4
        
        # restrict to simpler terms (tbm[5] = complicated diagonally-directed hopping term)
        tbm′ = tbm[1:4] 
        Δtbm′_m   = subduced_complement(tbm′, 10)                       # break mirror
        Δtbm′_tr  = subduced_complement(tbm′, 11; timereversal = false) # break TR
        Δtbm′_mtr = subduced_complement(tbm′, 10; timereversal = false) # break both
        @test length(Δtbm′_m) == 0
        @test length(Δtbm′_tr) == 0
        @test length(Δtbm′_mtr) == 1

        Δtbm′_C4 = subduced_complement(tbm′, 6)               # break C₄
        @test length(Δtbm′_C4) == 3
        Δtbm_diagonal_term = subduced_complement(tbm[5:5], 6) # break C₄ (diagonal term only)
        @test length(Δtbm_diagonal_term) == 1
        
        # adding more orbits, we find a simple mirror-breaking term, but not a TR-breaking
        tbm_big_full = tb_hamiltonian(cbr, [[0,0], [1,0], [1,1]])
        tbm_big = tbm_big_full[[1:4...,6]] # drop complicated diagonally-directed term again
        Δtbm_big_m   = subduced_complement(tbm_big, 10)                       # break mirror
        Δtbm_big_tr  = subduced_complement(tbm_big, 11; timereversal = false) # break TR
        Δtbm_big_mtr = subduced_complement(tbm_big, 10; timereversal = false) # break both
        @test length(Δtbm_big_m) == 1
        @test length(Δtbm_big_tr) == 0
        @test length(Δtbm_big_mtr) == 2
        @test issubset(Δtbm_big_m, Δtbm_big_mtr) # must subset eachother

        # breaking mirror and TR symmetry together should give the same basis as starting
        # directly with plane group p4 (#10) and breaking TR from the get-go
        brs10 = calc_bandreps(10, Val(D); timereversal=false)
        cbr10 = @composite brs10[1] # (2c|A) (unlike #11: two distinct M irreps, M₃ & M₄)

        tbm10 = tb_hamiltonian(cbr10, [[0,0], [1,0], [0,1]]) # (⋆)
        @test length(tbm10) == length(tbm) + length(Δtbm_mtr)
        # (⋆): must add [0,1] also, cf. diagonally-directed hopping term involving both 
        # [1,0] and [0,1] in C₄ setting

        tbm10_big = tb_hamiltonian(cbr10, [[0,0], [1,0], [1,1]])
        Δtbm_big_full_mtr = subduced_complement(tbm_big_full, 10; timereversal = false)
        @test length(tbm10_big) == length(tbm_big_full) + length(Δtbm_big_full_mtr)
    end

    @testset "3D example" begin
        # this is not a well-thought out example, but just added to test that things work
        # without erroring
        brs = calc_bandreps(16, Val(3); timereversal = true)
        cbr = @composite brs[1] + brs[end] #  (1h|A) + (1a|B₂) (2 bands)
        tbm = tb_hamiltonian(cbr, [[0,0,0],])

        @test length(subduced_complement(tbm, 3)) == 1
    end

    @testset "3D example, body-centered (I) lattice" begin
        brs = calc_bandreps(121, Val(3); timereversal = true)
        cbr = @composite brs[1] + brs[end-1] # (4d|A) + (2a|A₂) (3 bands)
        tbm = tb_hamiltonian(cbr, [[0,0,0],])

        Δtbm = subduced_complement(tbm, 82)
        @test length(Δtbm) == 2

        # completeness: `(4d|A)` of ⋕121 splits into `2c ⊕ 2d` in ⋕82, and `(2a|A₂)` maps
        # to one of `(2a|A/B)`; every such 3-band composite of ⋕82 has 6 = 4 + 2 terms
        brs82 = calc_bandreps(82, Val(3); timereversal = true)
        for (i, j, k) in Iterators.product((4, 5), (1, 2), (10, 11)) # 2c, 2d, 2a
            cbr82 = CompositeBandRep([n ∈ (i,j,k) ? 1 : 0 for n in eachindex(brs82)], brs82)
            @test length(tb_hamiltonian(cbr82, [[0,0,0],])) == length(tbm) + length(Δtbm)
        end

        # having added the complement, there is nothing further to find in ⋕82
        @test length(subduced_complement(vcat(tbm, Δtbm), 82)) == 0
    end

    @testset "composite band representations at a shared Wyckoff position" begin
        # terms belonging to *different* blocks can carry equal (`==`) hopping orbits, if
        # the associated band representations sit at the same Wyckoff position; such terms
        # must not be grouped together, since they do not share a coefficient basis
        brs = calc_bandreps(47, Val(3); timereversal = true) # P4/mmm
        cbr = @composite brs[57] + brs[60] # (1a|Ag) + (1a|B₁ᵤ)

        tbm = tb_hamiltonian(cbr, [[0,0,0]]) # on-site only: blocks (1,1) & (2,2), both δ=0
        @test length(tbm) == 2
        # subducing to G itself, with unchanged time-reversal, must give nothing new
        @test length(subduced_complement(tbm, 47)) == 0
        @test length(subduced_complement(tbm, 47; timereversal = false)) == 0

        Rs = [[0,0,0], [1,0,0], [0,1,0], [0,0,1]]
        tbm_nn = tb_hamiltonian(cbr, Rs)
        @test length(tbm_nn) == 9
        @test length(subduced_complement(tbm_nn, 47)) == 0
        Δtbm_nn = subduced_complement(tbm_nn, 47; timereversal = false)
        @test length(Δtbm_nn) == 1

        # the complement must be *complete*: breaking time-reversal in G should give the
        # same number of terms as building the model without time-reversal from the start
        brs′ = calc_bandreps(47, Val(3); timereversal = false)
        cbr′ = @composite brs′[57] + brs′[60]
        # ⋕57 & ⋕60 index the same band representations with and without time-reversal, but
        # that is a property of `calc_bandreps`' ordering rather than something we control
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
        tbm2 = tb_hamiltonian(cbr2, [[0,0,0], [1,0,0]])
        @test length(tbm2) == 6
        @test length(subduced_complement(tbm2, 47)) == 0
        @test length(subduced_complement(tbm2, 47; timereversal = false)) == 2
    end

    @testset "term grouping" begin
        # terms sharing a coefficient basis must be grouped together regardless of whether
        # they appear contiguously in the model (models may be built by `vcat` or indexing)
        brs2d = calc_bandreps(11, Val(2); timereversal = true)
        tbm = tb_hamiltonian((@composite brs2d[1]), [[0,0], [1,0]])
        @test _group_terms_by_block_and_orbit(tbm) == [[1], [2], [3, 4], [5]]

        p = [3, 1, 2, 4, 5] # splits the {3,4} group apart
        @test _group_terms_by_block_and_orbit(tbm[p]) == [[1, 4], [2], [3], [5]]
        for (sgnumᴴ, timereversal) in ((10, true), (6, true), (11, false), (10, false))
            @test length(subduced_complement(tbm[p], sgnumᴴ; timereversal)) ==
                  length(subduced_complement(tbm, sgnumᴴ; timereversal))
        end

        # equal-orbit terms from distinct blocks must stay in distinct groups, also when
        # they are adjacent
        brs = calc_bandreps(47, Val(3); timereversal = true)
        tbm2 = tb_hamiltonian((@composite brs[57] + brs[57]), [[0,0,0], [1,0,0]])
        q = [1, 3, 5, 2, 4, 6] # interleave, so that equal orbits become adjacent
        @test [t.block_ij for t in tbm2.terms] ==
              [(1,1), (1,1), (2,2), (2,2), (1,2), (1,2)]
        @test _group_terms_by_block_and_orbit(tbm2[q]) == [[i] for i in 1:6]
        @test length(subduced_complement(tbm2[q], 47; timereversal = false)) ==
              length(subduced_complement(tbm2,    47; timereversal = false))
    end

    @testset "centered lattices" begin
        # the subgroup generators must be converted to the primitive setting before the
        # constraints are imposed; if not, centered lattices error out in
        # `sgrep_induced_by_siteir` (which compares against primitivized site groups)
        brs = calc_bandreps(12, Val(3); timereversal = true) # C2/m (C-centered)
        cbr = @composite brs[1] # (4f|Ag)
        tbm = tb_hamiltonian(cbr, [[0,0,0], [1,0,0]])
        @test length(tbm) == 6

        # subducing to G itself, with unchanged time-reversal, must give nothing new
        @test length(subduced_complement(tbm, 12)) == 0

        @test length(subduced_complement(tbm, 12; timereversal = false)) == 1 # break TR
        @test length(subduced_complement(tbm, 5)) == 2                       # break mirror
        @test length(subduced_complement(tbm, 8)) == 2                       # break C₂ & -1
    end

    @testset "orbits that are empty in the parent group (issue #117)" begin
        # a (block, orbit) pair on which *every* coefficient is forbidden in G carries no
        # term, and so cannot be found without `Rs` - even though the symmetry reduction
        # may be exactly what allows it
        Rs = [[0,0,0], [1,0,0], [0,1,0], [0,0,1]]
        brs = calc_bandreps(47, Val(3); timereversal = true) # P4/mmm
        cbr = @composite brs[57] + brs[60] # (1a|Ag) + (1a|B₁ᵤ): s & p_z on a shared site
        tbm = tb_hamiltonian(cbr, Rs)
        @test length(tbm) == 9

        # the (1,2) block is forbidden on the x and y bonds in P4/mmm; dropping m_z (while
        # keeping inversion and TR) frees the bond along the retained 2-fold axis
        gensᴴ = [S"x,-y,-z", S"-x,-y,-z"] # 2ₓ & -1, i.e. 2/m with unique axis a
        @test length(_subduced_complement(tbm, gensᴴ)) == 0      # without `Rs`: misses it
        Δtbm = _subduced_complement(tbm, Rs, gensᴴ)              # over `Rs`: finds it
        @test length(Δtbm) == 1
        @test only(Δtbm).block_ij == (1, 2)
        @test length(subduced_complement(tbm, Rs, 10)) == 1      # ⋕10 = P2/m

        # completeness: the direct P2/m model over the same range has exactly one more term
        brs10 = calc_bandreps(10, Val(3); timereversal = true)
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

        # `Rs` must add the pairs that carry no term, and only those
        groups = _subduction_groups(tbm, Rs)
        @test sum(length(idxs) for (_, idxs) in groups) == length(tbm)
        @test count(((_, idxs),) -> isempty(idxs), groups) > 0
        @test [idxs for (_, idxs) in _subduction_groups(tbm, nothing)] ==
              _group_terms_by_block_and_orbit(tbm)
    end

    @testset "hopping range `Rs`" begin
        brs = calc_bandreps(11, Val(2); timereversal = true) # p4mm
        cbr = @composite brs[1] # (2c|A₁)
        Rs = [[0,0], [1,0]]
        tbm = tb_hamiltonian(cbr, Rs)

        # when the model already spans every orbit over `Rs`, both forms must agree
        for (sgnumᴴ, timereversal) in ((11, true), (10, true), (11, false), (10, false),
                                       (6, true), (6, false))
            @test length(subduced_complement(tbm, Rs, sgnumᴴ; timereversal)) ==
                  length(subduced_complement(tbm, sgnumᴴ; timereversal))
        end
        @test length(subduced_complement(tbm, Rs, 11)) == 0 # subducing to G itself

        # over `Rs`, a dropped term always comes back, whether or not the rest of its orbit
        # survived; without `Rs`, only the former (cf. issue #117)
        @test _group_terms_by_block_and_orbit(tbm) == [[1], [2], [3, 4], [5]]
        @test length(subduced_complement(tbm[[1,2,3,5]], Rs, 11)) == 1 # dropped term 4
        @test length(subduced_complement(tbm[[1,2,3,5]], 11)) == 1
        @test length(subduced_complement(tbm[1:4], Rs, 11)) == 1       # dropped term 5
        @test length(subduced_complement(tbm[1:4], 11)) == 0           # ← the asymmetry

        # a narrower `Rs` cannot report terms already in `tbm` as new
        @test length(subduced_complement(tbm, [[0,0]], 11)) == 0
        @test length(subduced_complement(tbm, [[0,0]], 10)) ==
              length(subduced_complement(tbm, 10))

        # a wider `Rs` also picks up longer-range terms
        Rs_big = [[0,0], [1,0], [1,1]]
        @test length(subduced_complement(tbm, Rs_big, 11)) ==
              length(tb_hamiltonian(cbr, Rs_big)) - length(tbm)

        # completeness against a directly-built subgroup model: the 2c orbit of p4mm splits
        # into 1c ⊕ 1b of p2mm. The orbits are those of p4mm, and reach further than `Rs`,
        # so the ⋕6 model must be built over a range covering the same hopping vectors
        brs6 = calc_bandreps(6, Val(2); timereversal = true)
        cbr6 = CompositeBandRep([n ∈ (5, 9) ? 1 : 0 for n in eachindex(brs6)], brs6)
        @test string(cbr6) == "(1c|A₁) + (1b|A₁)"
        @test length(tb_hamiltonian(cbr6, [[0,0], [1,0], [0,1], [-1,0]])) ==
              length(tbm) + length(subduced_complement(tbm, Rs, 6))

        # non-Hermitian models iterate over all blocks, not just the upper-triangular ones
        tbm_nh = tb_hamiltonian(cbr, Rs, Val(NONHERMITIAN))
        @test length(subduced_complement(tbm_nh, Rs, 11)) == 0
        for sgnumᴴ in (10, 6)
            @test length(subduced_complement(tbm_nh, Rs, sgnumᴴ)) ==
                  length(subduced_complement(tbm_nh, sgnumᴴ))
        end
    end

    @testset "subgroup precondition" begin
        # `_subduced_complement` asserts that `gensᴴ` generate a subgroup of G, given in G's
        # conventional setting; unreachable via `subduced_complement`, but check it can fire
        brs2d = calc_bandreps(11, Val(2); timereversal = true) # p4mm
        tbm = tb_hamiltonian((@composite brs2d[1]), [[0,0], [1,0]])

        @test _issubgroup([S"y,x"], 11)     # mₓᵧ ∈ p4mm
        @test !_issubgroup([S"-y,x+y"], 11) # 6⁺ ∉ p4mm
        @test length(_subduced_complement(tbm, [S"y,x"])) isa Int
        @test_throws AssertionError _subduced_complement(tbm, [S"-y,x+y"])
    end
end
