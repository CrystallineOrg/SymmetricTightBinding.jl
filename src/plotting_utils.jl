"""
    hop_plot(hop::HoppingOrbit{D}) where {D} --> Figure

Plots the hopping orbit `hop` in 2D or 3D. Also indicates the starting and ending atoms of the
hoppings and the bonds between them.
"""
function hop_plot(hop::HoppingOrbit{D}) where {D}
    P, V = Point{D,Float32}, Vec{D,Float32}
    starting_atoms = reduce(vcat, [[P(constant(h[1])) for h in hs] for hs in hop.hoppings])
    ending_atoms = reduce(vcat, [[P(constant(h[2] + h[3])) for h in hs] for hs in hop.hoppings])

    # figure / axis setup
    f = Figure()
    ax = if D == 3
        Axis3(f[1, 1]; aspect=:data, perspectiveness=0.2)
    else
        Axis(f[1, 1]; aspect=1)
    end
    hidespines!(ax)
    hidedecorations!(ax)

    # plot unit cell
    rect = Rect{D,Float32}(P(-0.5), V(1))  # unit cube at origin
    lines!(ax, rect; color=:black)

    # plot bonds
    arrows!(ax, starting_atoms, V.(ending_atoms - starting_atoms) .* 0.8;
        color=:gray,
        linewidth=0.025,
        arrowsize=D == 3 ? V(0.075, 0.075, 0.1) : V(0.075, 0.1), label="Bonds")

    # plot atoms
    meshscatter!(ax, unique(ending_atoms), markersize=0.05, color=:red, label="Ending atoms")
    meshscatter!(ax, unique(starting_atoms), markersize=0.05, color=:blue, label="Starting atoms")

    xmax = maximum(r -> abs(r[1]), ending_atoms)
    xlims!(ax, -xmax, xmax)
    ymax = maximum(r -> abs(r[2]), ending_atoms)
    ylims!(ax, -ymax, ymax)
    D == 3 && (zmax = maximum(r -> abs(r[3]), ending_atoms); zlims!(ax, -zmax, zmax))

    axislegend(ax)

    return f
end