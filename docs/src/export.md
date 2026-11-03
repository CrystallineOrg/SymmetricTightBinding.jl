# [Exporting to Wannier90](@id export)

[`write_wannier90`](@ref) writes a parameterized model in the file formats of [Wannier90](https://wannier.org), which many other tight-binding packages can read: e.g., [PythTB](https://www.physics.rutgers.edu/pythtb/), [TBmodels](https://tbmodels.greschd.ch) (and, from there, [Kwant](https://kwant-project.org)), [WannierTools](https://github.com/quanshengwu/wannier_tools), [WannierBerri](https://wannier-berri.org), and [WannierNLQG.jl](https://github.com/ZhuochengLu/WannierNLQG). The model is written via its real-space hopping matrices ``H(\mathbf{R})``, with elements ``H_{ij}(\mathbf{R}) = ⟨i,\mathbf{0}|\hat{H}|j,\mathbf{R}⟩``, to the [`_hr.dat` format](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/#seedname_hrdat).

## Exporting a simple model: graphene

We build a graphene model with nearest-neighbor hoppings:

```@example export
using Crystalline, SymmetricTightBinding

brs = bandreps(17, Val(2))
cbr = @composite brs[5] # (2b|A₁)
tbm = tb_hamiltonian(cbr)
ptbm = tbm([0.0, 1.0])
nothing # hide
```

To export the model, we pass a file prefix and the (primitive) lattice basis that the model refers to:

```@example export
lattice = directbasis(17, Val(2))
dir = mktempdir()
paths = write_wannier90(joinpath(dir, "graphene"), ptbm; lattice)
basename.(values(paths))
```

The `_hr.dat` file holds the hopping matrices ``H(\mathbf{R})``; the `_centres.xyz` and `.win` files hold the orbital positions (as Wannier centres) and the lattice vectors. As Wannier90's formats are three-dimensional, a two-dimensional model is placed in the ``xy``-plane:

```@example export
for (key, path) in pairs(paths)
    println("─── ", key, " ───")
    println(read(path, String))
end
```

The files can then be loaded e.g. in Python via PythTB or TBmodels:

```python
import pythtb
m_pythtb = pythtb.W90(dir, "graphene").model(zero_energy=0.0)

import tbmodels
m_tbmodels = tbmodels.Model.from_wannier_files(
    hr_file=f"{dir}/graphene_hr.dat",
    xyz_file=f"{dir}/graphene_centres.xyz",
    win_file=f"{dir}/graphene.win",
)
```
