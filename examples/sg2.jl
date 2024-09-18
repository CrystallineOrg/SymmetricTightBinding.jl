using Pkg
Pkg.activate(@__DIR__)

## necesary packages
using Crystalline
using SymmetryBases, MPBUtils
using PhotonicBandConnectivity
using TETB

sg_num = 2
brs = calc_bandreps(sg_num)
lgirsv = irreps(brs)

# ----------------------------------------------------------------------------------------#
s₁ = "[-Γ₁⁺+3Γ₁⁻,2R₁⁻, 2T₁⁻, 2U₁⁻, V₁⁺+V₁⁻, X₁⁺+X₁⁻,Y₁⁺+Y₁⁻,2Z₁⁺]"
v₁ = parse(SymmetryVector, s₁, lgirsv)

#s₂ =
# ----------------------------------------------------------------------------------------#

t = 1
d = stack(brs)[end, :]
long_modes = find_auxiliary_modes(t, d, brs)

### compute all possible decomposition into EBRs of vᵀ using the additional modes computed
all_band_repre = find_all_band_representations(v₁, long_modes, d, brs)