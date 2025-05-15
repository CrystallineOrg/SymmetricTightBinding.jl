using Pkg
Pkg.activate(@__DIR__)
using Crystalline, TETB

brs = calc_bandreps(17, Val(2));
coefs = zeros(Int, length(brs));
coefs[5] = 1;

cbr = CompositeBandRep(coefs, brs)

ptbm = tb_hamiltonian(cbr, [zeros(Int, dim(cbr))])([0.0, 1.0]);

using Brillouin, GLMakie

kpi = interpolate(irrfbz_path(17, directbasis(17, Val(2))), 100);

plot(kpi, spectrum(ptbm, kpi))
