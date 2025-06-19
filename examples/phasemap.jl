using Crystalline, SymmetricTightBinding 
using PythonCall
const pm = pyimport("phasemap") # requires `pkg> conda pip_add PhaseMap`

## --------------------------------------------------------------------------------------- #
# Define the tight-binding model we want to explore

brs = calc_bandreps(2, Val(3))
cbr = @composite brs[1] + brs[3]

tbm = tb_hamiltonian(cbr, [[0,0,0], [0,0,1]])
deleteat!(tbm.terms, [4])

## --------------------------------------------------------------------------------------- #
# Define some phase function over a two-dimensional parameter space

function phase(tbm, ϕθ)
    # set c1, c2, c3 to spherical coordinates
    ϕ, θ = ϕθ
    c1 = cos(ϕ) * cos(θ)
    c2 = cos(ϕ) * sin(θ)
    c3 = sin(ϕ)
    ptbm = tbm([0.25, c1, c2, c3, -0.25])
    collect_compatible(ptbm)
end
function phase_hash(tbm, ϕθ)
    hash(phase(tbm, ϕθ)) % Int
end
 

## --------------------------------------------------------------------------------------- #
# Run phasemap's algorithm to construct the phase diagram

limits = [float.((0, 2π)), float.((0, π))]
pm_phase_function = pyfunc(ϕθ->phase_hash(tbm, pyconvert(Tuple{Float64, Float64}, ϕθ)))
res = pm.run(pm_phase_function;
        limits=limits,
        mesh=6,
        num_steps=7)

## --------------------------------------------------------------------------------------- #
# We could now plot the results using phasemap's matplotlib functions, e.g., via this:
#       plt = pm.plot.plt # matplotlib.pyplot (equiv. to `pyimport("matplotlib.pyplot")`)
#       pm.plot.boxes(res)
#       plt.savefig("phasemap-example-sg2.pdf", bbox_inches="tight")
#
# But it's nice to have the results in Julia to have more control, so we extract box-data
# and then use GLMakie to visualize the phase diagram afterwards

# first, we extract all the info (box corners, sizes, phase values) from `res`
s_limits = (limits[1][1], limits[2][1]) # "origins" of limits
w_limits = (limits[1][2] - limits[1][1], limits[2][2] - limits[2][1]) # widths of limits
boxes_py = pyconvert(Vector, res.boxes)
box_corners = Vector{NTuple{2, Float64}}(undef, length(boxes_py))
box_sizes   = Vector{NTuple{2, Float64}}(undef, length(boxes_py))
box_phases  = Vector{Int}(undef, length(boxes_py))
for (i, bpy) in enumerate(boxes_py)
    scaled_corner = (pyconvert(Float64, bpy.corner[0]), pyconvert(Float64, bpy.corner[1])) # ∈[0,1]²
    scaled_size   = (pyconvert(Float64, bpy.size[0]), pyconvert(Float64, bpy.size[1]))     # ∈[0,1]²

    box_corners[i] = (s_limits[1] + w_limits[1]*scaled_corner[1], s_limits[2] + w_limits[2]*scaled_corner[2])
    box_sizes[i]  = (w_limits[1]*scaled_size[1], w_limits[2]*scaled_size[2])
    box_phases[i] = Bool(pytype(bpy.phase) == (@py int)) ? pyconvert(Int, bpy.phase) : 0
end

# then we change our arbitrary "phase values" from some hashed number to a consecutive
# number from 1 to N
u_phases_d = Dict(v=>i for (i,v) in enumerate(filter!(!iszero, unique(box_phases))))
remapped_box_phases = Vector{Float64}(undef, length(box_phases))
for (i, p) in enumerate(box_phases)
    remapped_box_phases[i] = iszero(p) ? NaN : u_phases_d[p]
end

# build a legend, using the symmetry vector of the valence band
legend_d = Dict{Int, String}()
for (i, (ϕθ, p)) in enumerate(zip(box_corners, box_phases))
    iszero(p) && continue
    id = u_phases_d[p]
    haskey(legend_d, id) && continue
    length(legend_d) == length(u_phases_d) && break
    n = phase(tbm, ϕθ)[1]
    s = replace(string(n), "Z₁⁻, Y₁⁻"=>"…", "T₁⁺, Γ₁⁺"=>"…", " (1 band)"=>"")
    t = calc_topology(n, brs)
    t == NONTRIVIAL && (s *= " (nontrivial)")
    legend_d[id] = s
end
legend = [legend_d[i] for i in 1:length(u_phases_d)]

## --------------------------------------------------------------------------------------- #
# Finally, we can plot the results using meshscatter with a sized box as marker

using GLMakie
f = Figure()
ax = Axis(f[1,1], 
    xticks = ([0, π, 2π], ["0", "π", "2π"]),   xlabel = "ϕ",
    yticks = ([0, π/2, π], ["0", "π/2", "π"]), ylabel = "θ",
    aspect = DataAspect()
    )
meshscatter!(ax, box_corners, markersize = box_sizes,
    marker = Rect2f(0, 0, 1, 1), 
    color = remapped_box_phases, 
    colormap= :tab10, colorrange = (1, 10),
    nan_color = :gray,
    label = [l=>(; color=i, markersize=.75) for (i,l) in enumerate(legend)])
xlims!(limits[1]...); ylims!(limits[2]...);
Legend(f[1,2], ax, framevisible=false)
f
