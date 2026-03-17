using Test
using Random
using SymmetricTightBinding
using Crystalline
using Crystalline: free

function _test_symmetry_analysis(brs, i, αβγ; rng_seed=1234)
    coefs = zeros(length(brs))
    coefs[i] = 1
    cbr = CompositeBandRep(coefs, brs)
    iszero(free(position(brs[i]))) || pin_free!(brs, i => αβγ)
    tbm = tb_hamiltonian(cbr)
    rng = seed!(copy(Random.default_rng()), rng_seed)
    ptbm = tbm(randn(rng, length(tbm)))
    ns = collect_compatible(ptbm)
    return !isempty(ns) && sum(ns) == SymmetryVector(cbr)
end

@testset "Symmetry analysis (full scan over EBRs)" begin
    # construct a tb-model for every EBR and verify that `collect_compatible` returns a set
    # of symmetry vectors whose sum is equal to the underlying EBR's: note that we cannot
    # generally require that `collect_compatible` returns a _single_ symmetry vector
    # equaling the EBR's symmetry vector, since the EBR might be decomposable (either into
    # fragile or genuinely topological EBRs): i.e., we don't know if the model is a single
    # connected band or not, only what the "sum" of symmetry vectors will be
    for D in 1:3
        αβγ = D == 1 ? [.1] : D == 2 ? [.1, .2] : [.1, .2, .3] # for `pin_free!`
        for sgnum in MAX_SGNUM[D]
            @testset "Space group $sgnum in dimension $D" begin
                brs = calc_bandreps(sgnum, Val(D))
                for i in eachindex(brs)
                    @test _test_symmetry_analysis(brs, i, αβγ)
                end
            end
        end
    end
end
