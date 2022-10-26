using GLMakie, Colors
include("./plant.jl")
include("./compute_curves.jl")
include("./scan.jl")

function cascade_plot_cashift!(ax, xs, gsyn, E, Iapp;
    res = 100, xbeg = 0, xend = 1, cabeg = 0, caend = 2)

    limits!(ax, cabeg, caend, xbeg, xend)
    l = ax.limits[]
    #define search space and scan
    Xsp = range(l[2][1], stop = l[2][2], length = res)
    Csp = range(l[1][1], stop = l[1][2], length = res)

    c1 = colorant"red"
    c2 = colorant"darkblue"
    N = 8 #number of curves
    colors = range(c1, stop = c2, length = N)

    v = collect(-100:.1:30)
    Ieq = get_fast_equilibrium(v, 0, 0, 0)

    for (i, cs) in enumerate(range(-100., stop = 350., length = N))
        # xshift, cashift, gsyn, E, Iapp. See definition for optional kw arguments.

        avg, freq, dxarr, dcarr = slowscan(Xsp,Csp,cs,xs, gsyn, E, Iapp)
        if i == 1
            #heatmap!(ax, Csp, Xsp, avg, colormap = :coolwarm,colorrange = (-70, -20))
            contour!(ax, Csp, Xsp, freq, levels = [0], color = :red, linewidth = 1)
            contour!(ax, Csp, Xsp, dxarr, levels = [0], color = :lightgrey,
            linewidth = 2)

            plot_x_nullcline_branches!(ax, 0, gsyn, E, Iapp; unstable_only = true, unstablecolor = :white, linewidth = 3)

        end
        if i == N
            _xs = -4
            avg, freq, dxarr, dcarr = slowscan(Xsp,Csp,cs,_xs, gsyn, E, Iapp)
            contour!(ax, Csp, Xsp, dxarr, levels = [0], color = :lightgrey,
            linewidth = 2)

            plot_x_nullcline_branches!(ax, -4, gsyn, E, Iapp; unstable_only = true, unstablecolor = :white, linewidth = 3)
        end
        contour!(ax, Csp, Xsp, dcarr, levels = [0], color = colors[i], linewidth = 3)
                #xnullc analytic
    end
    ax
end

"""F = let
     F = Figure()
     xs = 0; gsyn = 0; E = -70; Iapp = 0;
     ax = Axis(F[1,1])
     cascade_plot_cashift!(ax, xs, gsyn, E, Iapp; res = 80)
     F
end

save("cascade_plot_cashift.jpeg", F)"""
