import numpy as np
import pythtb as tb
import matplotlib.pyplot as plt

# Set up the lattice vectors
lat = [[1.0, 0.0], [-0.5, np.sqrt(3.0) / 2.0]]
# Set up the coordinates of the two atoms in the unit cell
orb = [[1.0 / 3.0, 2.0 / 3.0], [2.0 / 3.0, 1.0 / 3.0]]

# make a 2D tight-binding model
my_model = tb.tb_model(2, 2, lat, orb)

# set model parameters
delta = 0.0
t = -1.0

# set on-site energies and hoppings
my_model.set_onsite([-delta, delta])
my_model.set_hop(t, 0, 1, [0, 0])
my_model.set_hop(t, 0, 1, [-1, 0])
my_model.set_hop(t, 0, 1, [0, 1])

# visualize infinite model
(fig, ax) = my_model.visualize(0, 1)
ax.set_title("Graphene structure")
ax.set_xlabel("x coordinate")
ax.set_ylabel("y coordinate")
fig.tight_layout()
fig.savefig("graphene_structure.pdf")

# graphene band structure
path = [[0.0, 0.0], [1.0 / 3.0, 1.0 / 3.0], [1.0 / 2.0, 0.0], [0.0, 0.0]]
(k_vec, k_dist, k_node) = my_model.k_path(path, 100, report=False)
evals = my_model.solve_all(k_vec)

# plot band structure
fig, ax = plt.subplots()

for n in range(evals.shape[0]):
    ax.plot(k_dist, evals[n, :], color="black")

ax.set_title("Graphene band structure")
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Energy")
fig.tight_layout()
fig.savefig("graphene_band_structure.pdf")
