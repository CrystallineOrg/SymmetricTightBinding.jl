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

# ----------------------------------------------------------------------------------------#
# Symmetry vector test for different cases

s_1 = "[Γ₄⁻+Γ₁⁻+Γ₂⁺+Γ₅⁺, R₁⁻+R₂⁺+R₄⁻+R₅⁺ ,M₁+M₂+M₃+M₄,X₁+X₂+X₃+X₄]" # symmetry vector with
# possible non-physical solutions
v_1 = parse(SymmetryVector, s_1, lgirsv)

s_2 = "[-Γ₁⁺+2Γ₄⁻+Γ₂⁻, R₄⁻+R₅⁺,M₁+2M₄,X₁+X₃+X₄]" # "normal" symmetry vector
v_2 = parse(SymmetryVector, s_2, lgirsv)

s_3 = "[Γ₄⁻+Γ₂⁺+Γ₃⁻, R₂⁺+R₄⁻+R₃⁻+R₄⁻, M₁+M₂+M₄, X₁+X₂+X₃]" # needs high dim in auxiliary modes
v_3 = parse(SymmetryVector, s_3, lgirsv)
# ----------------------------------------------------------------------------------------#

t = 6
d = stack(brs)[end, :]
long_modes = find_auxiliary_modes(t, d, brs)

### compute all possible decomposition into EBRs of vᵀ using the additional modes computed
all_band_repre = find_all_band_representations(v_2, long_modes, d, brs)
