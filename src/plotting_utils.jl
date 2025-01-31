using GLMakie

function hopplot(hop::HoppingOrbit)
    starting_atoms = reduce(vcat, [[Point3f(constant(h[1])) for h in hs] for hs in hop.hoppings])
    ending_atoms = reduce(vcat, [[Point3f(constant(h[2] + h[3])) for h in hs] for hs in hop.hoppings])

    # figure / axis setup
    f = Figure()
    ax = Axis3(f[1, 1]; aspect=:data, perspectiveness=0.2,
        xspinesvisible=false, yspinesvisible=false, zspinesvisible=false,
        xticksvisible=false, yticksvisible=false, zticksvisible=false,
        xticklabelsvisible=false, yticklabelsvisible=false, zticklabelsvisible=false,
    )

    # plot unit cell
    rect = Rect3(Point3f(-0.5, -0.5, -0.5), Vec3f(1, 1, 1))  # Unit cube at (0,0,0)
    lines!(ax, rect; color=:black)

    # plot bonds
    for (a, b) in zip(starting_atoms, ending_atoms)
        #lines!(ax, [a, b], color = :gray, linewidth = 10)
        arrows!(ax, [a], [Vec(b - a)] * 0.8, color=:gray, linewidth=0.025,
            arrowsize=Vec3f(0.075, 0.075, 0.1))
    end

    # plot atoms
    end_atoms = meshscatter!(ax, unique(ending_atoms), markersize=0.05, color=:red,
        label="Ending atoms")
    start_atoms = meshscatter!(ax, unique(starting_atoms), markersize=0.05, color=:blue,
        label="Starting atoms")

    # indicate lengeds for atoms and bonds (need need to hardcode the markers -not implemented
    # on `Makie.jl` yet)
    axislegend(ax, framevisible=false)

    return f
end