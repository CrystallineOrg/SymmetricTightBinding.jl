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
s₁ = "[Γ₄⁻+Γ₁⁻+Γ₂⁺+Γ₅⁺, R₁⁻+R₂⁺+R₄⁻+R₅⁺ ,M₁+M₂+M₃+M₄,X₁+X₂+X₃+X₄]"
v₁ = parse(SymmetryVector, s₁, lgirsv)

#s₂ =
# ----------------------------------------------------------------------------------------#

t = 2
d = stack(brs)[end, :]
long_modes = find_auxiliary_modes(t, d, brs)

### compute all possible decomposition into EBRs of vᵀ using the additional modes computed
all_band_repre = find_all_band_representations(v₁, long_modes, d, brs)