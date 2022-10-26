using GLMakie

include("./plant.jl")
include("./compute_curves.jl")
include("./scan.jl")
include("./phase_plot.jl")
include("./cascade_plot_Iapp.jl")
include("./cascade_plot_cashift.jl")
include("./cascade_plot_xshift.jl")
include("./cascade_plot_gsyn.jl")

begin f = Figure(resolution=(1000,1500))

    # Row 1: sample fig, Iapp

    # replicating code from samplefig.jl

    res = 400;
    text_args = (textsize = 26, font = "TeX Gyre Heros Bold", padding = (0, 5, 5, 0), halign = :right)

    cs= 100; xs = -2; gsyn=0; E=-70; Iapp=0;
    ax2, sol1  = phase_plot!(Axis(f[2:6,1], ylabel="x-variable"), xs, cs, gsyn, E, Iapp; res = res, caend=2.5)
    hidexdecorations!(ax2, ticks=false)


    ax1 = Axis(f[1,1], title=L"Voltage \, Trace", xlabel="Time [sec]", ylabel="Voltage [mV]")
    ticks = collect(0:5000:10000)
    ticklabels = [ "$(round(x/1000))" for x in ticks ]
    ax1.xticks = (ticks, ticklabels)
    lines!(ax1, sol1.t, sol1[1,:], color = :black)

    ax = Axis(f[2:6,2], title=L"$I_{app}$",  ylabel="x-variable", titlesize = 22)
    cs = -60; xs = -4; gsyn = 0; E = -70;
    cascade_plot_Iapp!(ax, cs, xs, gsyn, E; res = res)
    hideydecorations!(ax, ticks=false)
    hidexdecorations!(ax, ticks=false)


    # Row 2: xshift, cashift
    cs = 0; gsyn = 0; E = -70; Iapp = 0;
    ax = Axis(f[7:11,1], title=L"$x$ shift",  ylabel="x-variable", titlesize = 22)
    cascade_plot_xshift!(ax, cs, gsyn, E, Iapp; res = res, caend=2.5, cabeg = 0)
    hidexdecorations!(ax, ticks=false)

    xs = 0; gsyn = 0; E = -70; Iapp = 0;
    ax = Axis(f[7:11,2], title=L"$Ca$ shift",  ylabel="x-variable", titlesize = 22)
    cascade_plot_cashift!(ax, xs, gsyn, E, Iapp; res = res, caend=2.5)
    hideydecorations!(ax, ticks=false)
    hidexdecorations!(ax, ticks=false)


    # Row 3: gsyn exc, gsyn inh

    cs = -60; xs = -4; gsyn = 0; E = -70; Iapp = 0;
    ax = Axis(f[12:16,1], title=L"$g_{syn}$ inhibition",xlabel="[Ca]-variable", ylabel="x-variable", titlesize = 22)
    cascade_plot_gsyn!(ax, cs, xs, E, Iapp; res = res)

    cs = -60; xs = -4; gsyn = 0; E = 30; Iapp = 0;
    ax = Axis(f[12:16,2], title=L"$g_{syn}$ excitation",xlabel="[Ca]-variable", ylabel="x-variable", titlesize = 22)
    cascade_plot_gsyn!(ax, cs, xs, E, Iapp; res = res, gmax = .002)
    hideydecorations!(ax, ticks=false)

    colgap!(f.layout, 20)
    rowgap!(f.layout, 1)

    f
end

save("cascade_plots.jpeg", f)
