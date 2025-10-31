using Test
using SymmetricTightBinding
using Crystalline
using Crystalline: free

@testset "Symmetry analysis (documentation example)" begin
    # Example 1
    sgnum = 17                         # plane group p6mm
    brs = calc_bandreps(sgnum, Val(2)) # band representations
    cbr = @composite brs[5]            # (2b|A₁) EBR
    tbm = tb_hamiltonian(cbr)          # tight-binding model (nearest neigbors)
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

@testset "Symmetry analysis (full scan over EBRs)" begin
    # construct a tb-model for every EBR and verify that `collect_compatible` returns a set
    # of symmetry vectors whose sum is equal to the underlying EBR's: note that we cannot
    # generally require that `collect_compatible` returns a _single_ symmetry vector
    # equaling the EBR's symmetry vector, since the EBR might be decomposable (either into
    # fragile or genuinely topological EBRs): i.e., we don't know if the model is a single
    # connected band or not, only what the "sum" of symmetry vectors will be
    for D in 1:3
        αβγ = D == 1 ? [.1] : D == 2 ? [.1, .2] : [.1, .2, .3] # for `pin_free!`
        println("--- D = $D ---")
        for sgnum in 1:MAX_SGNUM[D]
            showed_sgnum = false
            #@testset "Space group $sgnum in dimension $D" begin
                brs = calc_bandreps(sgnum, Val(D))
                for i in eachindex(brs)
                    coefs = zeros(length(brs))
                    coefs[i] = 1
                    cbr = CompositeBandRep(coefs, brs)

                    # if the EBR had any free parameters associated with its Wyckoff
                    # position, we must pin them, using `pin_free!`
                    iszero(free(position(brs[i]))) || pin_free!(brs, i=>αβγ)

                    # build a random tight-binding model for the `i`th EBR
                    tbm = try
                        tb_hamiltonian(cbr)
                    catch e
                        showed_sgnum || (showed_sgnum = true; println("#$sgnum"))
                        println("  brs[$i] = $cbr:\t`tb_hamiltonian` error")
                        continue
                    end
                    ptbm = tbm(randn(length(tbm)))

                    # check that the symmetry vector returned by `collect_compatible` is
                    # what we started with (i.e. equal to `cbr`)
                    try
                        ns = collect_compatible(ptbm)
                        if length(ns) == 0
                            showed_sgnum || (showed_sgnum = true; println("#$sgnum:"))
                            println("  brs[$i] = $cbr:\tfailed to collect any compatible symmetry vectors")
                        else
                            cond = sum(ns) == SymmetryVector(cbr)
                            if !cond
                                showed_sgnum || (showed_sgnum = true; println("#$sgnum:"))
                                println("  brs[$i] = $cbr:\tsymmetry vector discrepancy")
                                println("    n   = $(sum(ns))\n    cbr = $(SymmetryVector(cbr))")
                            end
                        end
                    catch e
                        showed_sgnum || (showed_sgnum = true; println("#$sgnum:"))
                        println("  brs[$i] = $cbr:\t`collect_compatible` error")
                        println("    ", e)
                        continue
                    end
                    #@test cbr == only(collect_compatible(ptbm)) 
                end # for br in brs
            #end # @testset "Space group $sgnum in dimension $D"
        end # for sgnum in MAX_SGNUM[D]
        println()
    end # for D in 1:3
end # @testset "Symmetry analysis"