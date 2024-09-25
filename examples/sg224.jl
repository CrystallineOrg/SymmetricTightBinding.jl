using Pkg
Pkg.activate(@__DIR__)

### necessary packages
using Crystalline
using SymmetryBases, MPBUtils
using PhotonicBandConnectivity
using TETB

sg_num = 224
brs = calc_bandreps(sg_num)
lgirsv = irreps(brs)

#------------------------------------------------------------------------------------------#
# First minimal solutions in SG #224

s1 = "[Γ₄⁻+Γ₂⁻, R₂⁻+R₄⁻, M₁+M₄, X₁+X₃]"
m1 = parse(SymmetryVector, s1, lgirsv)

s2 = "[Γ₄⁻+Γ₂⁻, R₁⁺+R₅⁺, M₁+M₄, X₁+X₄]"
m2 = parse(SymmetryVector, s2, lgirsv)

s3 = "[-Γ₁⁺+Γ₄⁻+Γ₁⁻+Γ₂⁻, R₂⁺+R₄⁺, M₂+M₄, X₂+X₃]"
m3 = parse(SymmetryVector, s3, lgirsv)

s4 = "[-Γ₁⁺+Γ₄⁻+Γ₁⁻+Γ₂⁻, R₁⁻+R₅⁻, M₂+M₄, X₂+X₄]"
m4 = parse(SymmetryVector, s4, lgirsv)

#------------------------------------------------------------------------------------------#

μᴸ = 4
idxsᴸs = find_auxiliary_modes(μᴸ, brs)

### compute all possible decomposition into EBRs of m using the additional modes computed
candidatesv = find_apolar_modes(m3, idxsᴸs, brs)