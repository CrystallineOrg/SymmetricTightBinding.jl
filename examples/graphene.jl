using Pkg
Pkg.activate(@__DIR__)
using Crystalline, SymmetricTightBinding
using Brillouin, GLMakie

sgnum, D = 17, 2
brs = calc_bandreps(sgnum, Val(D))
cbr = @composite brs[5]

tbm = tb_hamiltonian(cbr)

ptbm = tbm(randn(length(tbm)))

Rs = directbasis(sgnum, Val(2))

kp = irrfbz_path(sgnum, Rs)
kpi = interpolate(kp, 200) # aim for 200 interpolations points

Es = spectrum(ptbm, kpi); # a 200×2 Matrix

plot(kpi, Es; annotations = collect_irrep_annotations(ptbm))