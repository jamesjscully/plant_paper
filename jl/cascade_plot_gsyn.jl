
using GLMakie, Colors
include("./plant.jl")
include("./compute_curves.jl")
include("./scan.jl")

function cascade_plot_gsyn!(ax, cs, xs, E, Iapp;
    res = 100, xbeg = 0, xend = 1, cabeg = 0, caend = 2.5, gmax = .01)

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
    for (i, gsyn) in enumerate(range(0., stop = gmax, length = N))
        # xshift, cashift, gsyn, E, Iapp. See definition for optional kw arguments.

        avg, freq, dxarr, dcarr = slowscan(Xsp,Csp,cs,xs, gsyn, E, Iapp)
        contour!(ax, Csp, Xsp, freq, levels = [0f0], color = colors[i], linewidth = 1.5)
        contour!(ax, Csp, Xsp, dxarr, levels = [0f0], color = colors[i], linewidth = 3)
        plot_x_nullcline_branches!(ax, xs, gsyn, E, Iapp;
            stablecolor = colors2[i], unstablecolor =colors2[i], linewidth = 3, unstable_only = true)
    end
    ax
end

"""F = let
    F = Figure()
    cs = -60; xs = -4; gsyn = 0; E = -70; Iapp = 0;
    ax = Axis(F[1,1])
    cascade_plot_gsyn!(ax, cs, xs, E, Iapp; res = 80)

    F
end
save("cascade_plot_gsyn_inh.jpeg", F)

F2 = let
    F = Figure()
    cs = -60; xs = -4; gsyn = 0; E = 30; Iapp = 0;
    ax = Axis(F[1,1])
    cascade_plot_gsyn!(ax, cs, xs, E, Iapp; res = 80, gmax = .001)
    F
end
save("cascade_plot_gsyn_exc.jpeg", F2)"""
