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
# First minimal solutions in SG #2

s_1 = "[-Γ₁⁺+3Γ₁⁻,2R₁⁻, 2T₁⁻, 2U₁⁻, V₁⁺+V₁⁻, X₁⁺+X₁⁻,Y₁⁺+Y₁⁻,2Z₁⁺]"
v_1 = parse(SymmetryVector, s_1, lgirsv)

s_2 = "[-Γ₁⁺+3Γ₁⁻, 2R₁⁻, T₁⁺ + T₁⁻, U₁⁺ + U₁⁻, V₁⁺ + V₁⁻, 2X₁⁻, 2Y₁⁻, 2Z₁⁺]"
v_2 = parse(SymmetryVector, s_2, lgirsv)

s_3 = "[-Γ₁⁺+3Γ₁⁻, R₁⁺+R₁⁻, T₁⁺ + T₁⁻, 2 U₁⁻, 2 V₁⁻, X₁⁺+X₁⁻, 2Y₁⁻, 2Z₁⁺]"
v_3 = parse(SymmetryVector, s_3, lgirsv)

s_4 = "[-Γ₁⁺+3Γ₁⁻, R₁⁺+R₁⁻, 2 T₁⁻, U₁⁺ + U₁⁻, 2 V₁⁻, 2X₁⁻, Y₁⁺+Y₁⁻, 2Z₁⁺]"
v_4 = parse(SymmetryVector, s_4, lgirsv)

s_5 = "[-Γ₁⁺+3Γ₁⁻, 2R₁⁻, T₁⁺ + T₁⁻, 2 U₁⁻, 2 V₁⁺, 2X₁⁻, Y₁⁺+Y₁⁻, Z₁⁺+Z₁⁻]"
v_5 = parse(SymmetryVector, s_5, lgirsv)

# ----------------------------------------------------------------------------------------#

t = 1
d = stack(brs)[end, :]
long_modes = find_auxiliary_modes(t, d, brs)

### compute all possible decomposition into EBRs of vᵀ using the additional modes computed
all_band_repre = find_all_band_representations(v_5, long_modes, d, brs)