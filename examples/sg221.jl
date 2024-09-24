using Pkg
Pkg.activate(@__DIR__)

### necessary packages
using Crystalline
using SymmetryBases, MPBUtils
using PhotonicBandConnectivity
using TETB

sg_num = 221
brs = calc_bandreps(sg_num)
lgirsv = irreps(brs)

#------------------------------------------------------------------------------------------#
# First minimal solutions in SG #221

s1 = "[-Γ₁⁺+Γ₄⁻, R₃⁻, M₂⁻+M₃⁺, X₅⁺]"
m1 = parse(SymmetryVector, s1, lgirsv)

s2 = "[-Γ₁⁺+Γ₄⁻, R₃⁺, M₂⁺+M₃⁻, X₅⁻]"
m2 = parse(SymmetryVector, s2, lgirsv)

#------------------------------------------------------------------------------------------#

μᴸ = 1
idxsᴸs = find_auxiliary_modes(μᴸ, brs)

### compute all possible decomposition into EBRs of m using the additional modes computed
candidates = find_apolar_modes(m2, idxsᴸs, brs)