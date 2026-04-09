using Test
using SymmetricTightBinding
using Crystalline
using DeepDiffs: deepdiff
using SymmetricTightBinding: TightBindingElementString

# test print with nicely printed diff on failures
# (adapted from Crystalline.jl's test/show.jl)
function test_show(expected::AbstractString, observed::AbstractString)
    if expected == observed
        @test true
    else
        old = Base.have_color
        @eval Base have_color = true
        try
            println(deepdiff(expected, observed))
        finally
            @eval Base have_color = $old
        end
        @test :expected == :observed
    end
end
test_tp_show(v, expected::AbstractString) = test_show(repr(MIME"text/plain"(), v), expected)

@testset "Show methods" begin
    # set up graphene model (plane group 17, (2b|A₁) EBR)
    brs = calc_bandreps(17, Val(2))
    br = brs[5]
    cbr = @composite brs[5]
    hop_orbits = obtain_symmetry_related_hoppings([[0, 0]], br, br)
    tbm = tb_hamiltonian(cbr, [[0, 0]])
    ptbm = tbm([0.0, 1.0])

    @testset "HoppingOrbit" begin
        str = """
        HoppingOrbit{2} (a + δ = b + R):
         δ₁ = [0, 0]: [([1/3, 2/3] → [1/3, 2/3] + [0, 0]), ([2/3, 1/3] → [2/3, 1/3] + [0, 0])]"""
        test_tp_show(hop_orbits[1], str)

        str = """
        HoppingOrbit{2} (a + δ = b + R):
         δ₁ = [-1/3, 1/3]:  [([2/3, 1/3] → [1/3, 2/3] + [0, 0])]
         δ₂ = [-1/3, -2/3]: [([2/3, 1/3] → [1/3, 2/3] + [0, -1])]
         δ₃ = [2/3, 1/3]:   [([2/3, 1/3] → [1/3, 2/3] + [1, 0])]
         δ₄ = [1/3, -1/3]:  [([1/3, 2/3] → [2/3, 1/3] + [0, 0])]
         δ₅ = [1/3, 2/3]:   [([1/3, 2/3] → [2/3, 1/3] + [0, 1])]
         δ₆ = [-2/3, -1/3]: [([1/3, 2/3] → [2/3, 1/3] + [-1, 0])]"""
        test_tp_show(hop_orbits[2], str)
    end

    @testset "TightBindingTerm" begin
        str = """
        2×2 TightBindingTerm{2} (hermitian) over [(2b|A₁)]:
         1  0
         0  1"""
        test_tp_show(tbm[1], str)

        str = """
        2×2 TightBindingTerm{2} (hermitian) over [(2b|A₁)]:
         0                     𝕖(-δ₄)+𝕖(-δ₅)+𝕖(-δ₆)
         𝕖(-δ₁)+𝕖(-δ₂)+𝕖(-δ₃)  0                   
        δ₁=[-1/3,1/3], δ₂=[-1/3,-2/3], δ₃=[2/3,1/3], δ₄=-δ₁, δ₅=-δ₂, δ₆=-δ₃"""
        test_tp_show(tbm[2], str)
    end

    @testset "TightBindingModel" begin
        str = """
        2-term 2×2 TightBindingModel{2} (hermitian) over (2b|A₁):
        ┌─
        1. ⎡ 1  0 ⎤
        │  ⎣ 0  1 ⎦
        └─ (2b|A₁) self-term
        ┌─
        2. ⎡ 0                     𝕖(-δ₄)+𝕖(-δ₅)+𝕖(-δ₆) ⎤
        │  ⎣ 𝕖(-δ₁)+𝕖(-δ₂)+𝕖(-δ₃)  0                    ⎦
        └─ (2b|A₁) self-term:  δ₁=[-1/3,1/3], δ₂=[-1/3,-2/3], δ₃=[2/3,1/3], δ₄=-δ₁, δ₅=-δ₂, δ₆=-δ₃"""
        test_tp_show(tbm, str)
    end

    @testset "ParameterizedTightBindingModel" begin
        str = """
        2-term 2×2 ParameterizedTightBindingModel{2} (hermitian) over (2b|A₁) with amplitudes:
         [0, 1.0]"""
        test_tp_show(ptbm, str)
    end

    @testset "TightBindingElementString" begin
        context = :color => true # enable color (ANSI codes) in `sprint` below
        @test sprint(show, TightBindingElementString("c₁e(δ₁)", true);  context) == "\e[34mc₁e(δ₁)\e[39m" # `active = true`  => blue
        @test sprint(show, TightBindingElementString("c₁e(δ₁)", false); context) == "\e[0mc₁e(δ₁)"        # `active = false` => white/default
        @test sprint(show, TightBindingElementString("0",       false); context) == "\e[90m0\e[39m"       # `active = false` && zero-string => light_dark
        @test sprint(show, TightBindingElementString("0",       true);  context) == "\e[0m0"              # `active = true`  && zero-string => white/default
    end

    @testset "summary" begin
        test_show(sprint(summary, tbm),
                  "2-term 2×2 TightBindingModel{2} (hermitian) over (2b|A₁)")
        test_show(sprint(summary, ptbm),
                  "2-term 2×2 ParameterizedTightBindingModel{2} (hermitian) over (2b|A₁)")
        test_show(sprint(summary, tbm[1]),
                  "2×2 TightBindingTerm{2} (hermitian) over [(2b|A₁)]")
    end
end
