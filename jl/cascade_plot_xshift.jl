using GLMakie, Colors
include("./plant.jl")
include("./compute_curves.jl")
include("./scan.jl")


function cascade_plot_xshift!(ax, cs, gsyn, E, Iapp;
    res = 100, xbeg = 0, xend = 1, cabeg = .4, caend = 1.5)

    limits!(ax, cabeg, caend, xbeg, xend)
    l = ax.limits[]
    #define search space and scan
    Xsp = range(l[2][1], stop = l[2][2], length = res)
    Csp = range(l[1][1], stop = l[1][2], length = res)

    c1 = colorant"red"
    c2 = colorant"darkblue"
    c3 = colorant"rgb(255,200,200)"
    c4 = colorant"lightblue"
    N = 8 #number of curves
    colors = range(c1, stop = c2, length = N)
    colors2 = range(c3, stop = c4, length = N)

    for (i, xshift) in enumerate(range(0., stop = -4., length = N))
        # xshift, cashift, gsyn, E, Iapp. See definition for optional kw arguments.

        avg, freq, dxarr, dcarr = slowscan(Xsp,Csp,cs,xshift, gsyn, E, Iapp)
        if i == 1
            #heatmap!(ax, Csp, Xsp, avg, colormap = :coolwarm,colorrange = (-70, -20))
            contour!(ax, Csp, Xsp, freq, levels = [0f0], color = :red, linewidth = 1)
        end
        contour!(ax, Csp, Xsp, dxarr, levels = [0f0], color = colors[i], linewidth = 3)
        plot_x_nullcline_branches!(ax, xshift, gsyn, E, Iapp;
            stablecolor = colors2[i], unstablecolor =colors2[i], unstable_only = true, linewidth = 3)
    end
    ax
end

"""F = let
    F = Figure()
    cs = 0; gsyn = 0; E = -70; Iapp = 0;
    ax = Axis(F[1,1])
    cascade_plot_xshift!(ax, cs, gsyn, E, Iapp; res = 80)
    F
end
F
save("cascade_plot_xshift.jpeg", F)"""
