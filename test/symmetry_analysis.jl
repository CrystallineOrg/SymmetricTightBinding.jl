using Test
using Random
using SymmetricTightBinding
using Crystalline
using Crystalline: free

# pre-existing failures in centered lattices (see devdocs/symmetry_eigenvalue_conventions.md)
const KNOWN_FAILURES = Dict(
    # SG => Set of broken EBR indices
    68  => Set(1:14),           # Ccca: all EBRs fail (C-centered)
    88  => Set([5, 6, 8, 9]),   # I4₁/a: 4a and 4b positions
    141 => Set([9,10,11,12,14,15,16,17]), # I4₁/amd: 4a and 4b (excl. E irreps)
    142 => Set([1, 2, 5, 8]),   # I4₁/acd: 16e and 8b positions
    214 => Set([1,4,5,8,9,10,12,13]), # I4₁32: 12c, 12d, 8a, 8b positions
    220 => Set([3, 4, 6, 7]),   # I-43d: 12a and 12b positions
    230 => Set([4, 7, 8, 9]),   # Ia-3d: 24c and 16b positions
)

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

max_sgnum_3d = 230 # NB: set low if exploring failures: "high" `sgnum`s take a lot of time
@testset "Symmetry analysis (full scan over EBRs)" begin
    # construct a tb-model for every EBR and verify that `collect_compatible` returns a set
    # of symmetry vectors whose sum is equal to the underlying EBR's: note that we cannot
    # generally require that `collect_compatible` returns a _single_ symmetry vector
    # equaling the EBR's symmetry vector, since the EBR might be decomposable (either into
    # fragile or genuinely topological EBRs): i.e., we don't know if the model is a single
    # connected band or not, only what the "sum" of symmetry vectors will be
    for D in 1:3
        αβγ = D == 1 ? [.1] : D == 2 ? [.1, .2] : [.1, .2, .3] # for `pin_free!`
        for sgnum in 1:min(MAX_SGNUM[D], max_sgnum_3d)
            @testset "Space group $sgnum in dimension $D" begin
                brs = calc_bandreps(sgnum, Val(D))
                for i in eachindex(brs)
                    is_known_failure = i in get(KNOWN_FAILURES, sgnum, Set{Int}())
                    if is_known_failure
                        @test_broken _test_symmetry_analysis(brs, i, αβγ)
                    else
                        @test _test_symmetry_analysis(brs, i, αβγ)
                    end
                end
            end
        end
    end
end
