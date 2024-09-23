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

s_1 = "[-Γ₁⁺+Γ₄⁻, R₃⁻, M₂⁻+M₃⁺, X₅⁺]"
v_1 = parse(SymmetryVector, s_1, lgirsv)

s_2 = "[-Γ₁⁺+Γ₄⁻, R₃⁺, M₂⁺+M₃⁻, X₅⁻]"
v_2 = parse(SymmetryVector, s_2, lgirsv)

#------------------------------------------------------------------------------------------#

t = 1
d = stack(brs)[end, :]
long_modes = find_auxiliary_modes(t, d, brs)

### compute all possible decomposition into EBRs of vᵀ using the additional modes computed
all_band_repre = find_all_band_representations(v_2, long_modes, d, brs)