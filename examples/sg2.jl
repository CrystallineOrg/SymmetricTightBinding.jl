using Pkg
Pkg.activate(@__DIR__)

# necessary packages
using Crystalline
using TETB

sgnum = 2
brs = calc_bandreps(sgnum)
lgirsv = irreps(brs)

# ----------------------------------------------------------------------------------------#
# First minimal solutions in SG #2

s1 = "[-Γ₁⁺+3Γ₁⁻,2R₁⁻, 2T₁⁻, 2U₁⁻, V₁⁺+V₁⁻, X₁⁺+X₁⁻,Y₁⁺+Y₁⁻,2Z₁⁺]"
m1 = parse(SymmetryVector, s1, lgirsv)

s2 = "[-Γ₁⁺+3Γ₁⁻, 2R₁⁻, T₁⁺ + T₁⁻, U₁⁺ + U₁⁻, V₁⁺ + V₁⁻, 2X₁⁻, 2Y₁⁻, 2Z₁⁺]"
m2 = parse(SymmetryVector, s2, lgirsv)

s3 = "[-Γ₁⁺+3Γ₁⁻, R₁⁺+R₁⁻, T₁⁺ + T₁⁻, 2 U₁⁻, 2 V₁⁻, X₁⁺+X₁⁻, 2Y₁⁻, 2Z₁⁺]"
m3 = parse(SymmetryVector, s3, lgirsv)

s4 = "[-Γ₁⁺+3Γ₁⁻, R₁⁺+R₁⁻, 2 T₁⁻, U₁⁺ + U₁⁻, 2 V₁⁻, 2X₁⁻, Y₁⁺+Y₁⁻, 2Z₁⁺]"
m4 = parse(SymmetryVector, s4, lgirsv)

s5 = "[-Γ₁⁺+3Γ₁⁻, 2R₁⁻, T₁⁺ + T₁⁻, 2 U₁⁻, 2 V₁⁺, 2X₁⁻, Y₁⁺+Y₁⁻, Z₁⁺+Z₁⁻]"
m5 = parse(SymmetryVector, s5, lgirsv)

# ----------------------------------------------------------------------------------------#

μᴸ = 1
idxsᴸs = find_auxiliary_modes(μᴸ, brs)

# compute all possible decomposition into EBRs of m using the additional modes computed
candidatesv = find_apolar_modes(m5, idxsᴸs, brs)