using Test
using Random
using SymmetricTightBinding
using Crystalline
using Crystalline: free

# ---------------------------------------------------------------------------------------- #

@testset "Symmetry analysis (documentation example)" begin
    # Example 1
    sgnum = 17                         # plane group p6mm
    brs = calc_bandreps(sgnum, Val(2)) # band representations
    cbr = @composite brs[5]            # (2b|A₁) EBR
    tbm = tb_hamiltonian(cbr)          # tight-binding model (nearest neighbors)
    ptbm = tbm([0, 1])                 # zero self-energy, nonzero nearest-neighbor hopping
    ns = collect_compatible(ptbm)
    @test only(ns) == SymmetryVector(cbr) == SymmetryVector(CompositeBandRep(ptbm))

    # Example 2
    cbr′ = @composite brs[3] + brs[5] # (2a|A₁) + (3c|B₂)
    tbm′ = tb_hamiltonian(cbr′)
    ptbm′ = tbm′([2.5, 0, 0.2, 0, -1, 0])
    ns′ = collect_compatible(ptbm′)
    @test length(ns′) == 1 # bands are overlapping and not separable
    @test only(ns′) == SymmetryVector(cbr′)

    ptbm′′ = tbm′([2.5, 0, 0.2, 0, -1, .5]) # turn on hybridization and split bands
    ns′′ = collect_compatible(ptbm′′)
    @test length(ns′′) == 2 # bands are overlapping and not separable
    @test ns′′[1] ∉ brs[[3, 5]]
    @test ns′′[2] ∉ brs[[3, 5]]
    @test sum(ns′′) == SymmetryVector(cbr′)
end

# ---------------------------------------------------------------------------------------- #

# Construct a tb-model for every EBR and verify that `collect_compatible` returns a set of
# symmetry vectors whose sum is equal to the underlying EBR's: note that we cannot generally
# require that `collect_compatible` returns a _single_ symmetry vector equaling the EBR's
# symmetry vector, since the EBR might be decomposable (either into fragile or genuinely
# topological EBRs): i.e., we don't know if the model is a single connected band or not,
# only what the "sum" of symmetry vectors will be.

function _test_symmetry_analysis(brs, i, αβγ; rng = seed!(copy(Random.default_rng()), 1234))
    coefs = zeros(length(brs))
    coefs[i] = 1
    cbr = CompositeBandRep(coefs, brs)
    iszero(free(position(brs[i]))) || pin_free!(brs, i => αβγ)
    tbm = tb_hamiltonian(cbr)
    ptbm = tbm(randn(rng, length(tbm)))
    ns = collect_compatible(ptbm)
    @test !isempty(ns)
    @test sum(ns) == SymmetryVector(cbr)
end

@testset "Symmetry analysis (full scan over EBRs)" begin
    for D in 1:3
        αβγ = D == 1 ? [.1] : D == 2 ? [.1, .2] : [.1, .2, .3] # for `pin_free!`
        for sgnum in MAX_SGNUM[D]
            @testset "Space group $sgnum in dimension $D" begin
                brs = calc_bandreps(sgnum, Val(D))
                for i in eachindex(brs)
                    _test_symmetry_analysis(brs, i, αβγ)
                end
            end
        end
    end
end
