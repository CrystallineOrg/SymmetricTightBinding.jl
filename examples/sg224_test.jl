using Pkg
Pkg.activate(@__DIR__)

## necesary packages
using Crystalline
using SymmetryBases, MPBUtils
using PhotonicBandConnectivity
using TETB

sg_num = 224
brs = calc_bandreps(sg_num)
lgirsv = irreps(brs)

## ---------------------------------------------------------------------------------------#
# Symmetry vector test for different cases

# symmetry vector w/ possible non-physical solutions
s1 = "[Γ₄⁻+Γ₁⁻+Γ₂⁺+Γ₅⁺, R₁⁻+R₂⁺+R₄⁻+R₅⁺ ,M₁+M₂+M₃+M₄,X₁+X₂+X₃+X₄]"
m1 = parse(SymmetryVector, s1, lgirsv)

# "normal" symmetry vector
s2 = "[-Γ₁⁺+2Γ₄⁻+Γ₂⁻, R₄⁻+R₅⁺,M₁+2M₄,X₁+X₃+X₄]"
m2 = parse(SymmetryVector, s2, lgirsv)

# needs high dim in auxiliary modes
s3 = "[Γ₄⁻+Γ₂⁺+Γ₃⁻, R₂⁺+R₄⁻+R₃⁻+R₄⁻, M₁+M₂+M₄, X₁+X₂+X₃]"
m3 = parse(SymmetryVector, s3, lgirsv)

# a mode that isn't connected to zero frequency; just a plain EBR (4c|Eᵤ) [must run with
# `connected_to_zero_frequency=false`]
s4 = "[M₁+M₂+M₃+M₄, X₁+X₂+X₃+X₄, Γ₃⁻+Γ₄⁻+Γ₅⁻, R₃⁺+R₄⁺+R₅⁺]" # (4c|Eᵤ) = brs[15] # BROKEN - see TODO in `find_bandrep_decompositions`
m4 = parse(SymmetryVector, s4, lgirsv)

## ----------------------------------------------------------------------------------------#

μᴸ = 2
idxsᴸs = find_auxiliary_modes(μᴸ, brs)

### compute all possible decomposition into EBRs of m using the additional modes computed
candidatesv = find_apolar_modes(m4, idxsᴸs, brs)

candidatesv′ = find_bandrep_decompositions(m3, brs; μᴸ_min=μᴸ)
