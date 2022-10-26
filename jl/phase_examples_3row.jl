using GLMakie
include("./plant.jl")
include("./compute_curves.jl")
include("./scan.jl")
include("./phase_plot.jl")

function formatphase!(ax1, ax2, args; res = 80, tmax = 20000, kwargs...)
  ax1, sol  = phase_plot!(ax1,args...; res = res, tmax = tmax, kwargs...)
  hidespines!(ax1)
  #hidedecorations!(ax1)
  lines!(ax2, sol.t, sol[1,:], color = :black)
  hidespines!(ax2)
  hidedecorations!(ax2, label = false)

  ax1, ax2
end


begin F = Figure(resolution=(1500,1500))
# xshift, cashift, gsyn, E, Iapp. See definition for optional kw arguments.
  res = 300
  text_args = (textsize = 26, font = "TeX Gyre Heros Bold", padding = (0, 0, 0, 0), halign = :right)
  axargs = (xticks = 2, yticks = 2)

  args = (xs = 0, cs = -60, gsyn = 0, E = -70, Iapp = 0)
  kwargs = (cabeg = .4, caend = .95)
  formatphase!(Axis(F[2:6,1]), Axis(F[1,1]),args; res = res, tmax = 20000,
    ca0 = 1.2, x0 = .6, kwargs...)
  rowgap!(F.layout,1, Relative(.01))

  args = (xs = 0, cs = 0, gsyn = 0, E = -70, Iapp = 0)
  formatphase!(Axis(F[2:6,2]), Axis(F[1,2]), args; res = res, tmax = 20000,
    ca0 = 1.2, x0 = .6, kwargs...)

  args = (xs = 0, cs = 120, gsyn = 0, E = -70, Iapp = 0)
  formatphase!(Axis(F[2:6,3]), Axis(F[1,3]), args; res = res, tmax = 20000,
    ca0 = .41, x0 = .3, kwargs...)
##
    kwargs = (cabeg = .4, caend = 1.1)
    args = (xs = -1.4, cs = -60, gsyn = 0, E = -70, Iapp = 0)
    formatphase!(Axis(F[8:12,1]), Axis(F[7,1]), args; res = res, tmax = 20000,
      ca0 = .55, x0 = .1, kwargs...)
    rowgap!(F.layout,7, Relative(.01))


    args = (xs = -1.4, cs = -32, gsyn = 0, E = -70, Iapp = 0)
    formatphase!(Axis(F[8:12,2]), Axis(F[7,2]), args; res = res, tmax = 60000,
      ca0 = .55, x0 = .1, kwargs...)

    args = (xs = -1.4, cs =  40, gsyn = 0, E = -70, Iapp = 0)
    formatphase!(Axis(F[8:12,3]), Axis(F[7,3]), args; res = res, tmax = 20000,
      ca0 = .55, x0 = .1, kwargs...)
##
  kwargs = (cabeg = .4, caend = 1.5)
  args = (xs = -4, cs = -60, gsyn = 0, E = -70, Iapp = 0)
  formatphase!(Axis(F[14:18,1]), Axis(F[13,1]), args; res = res, tmax = 20000,
    ca0 = .55, x0 = .1, kwargs...)
  rowgap!(F.layout,7, Relative(.01))


  args = (xs = -4, cs = 0, gsyn = 0, E = -70, Iapp = 0)
  formatphase!(Axis(F[14:18,2]), Axis(F[13,2]), args; res = res, tmax = 20000,
    ca0 = .55, x0 = .1, kwargs...)

  args = (xs = -4, cs = 125, gsyn = 0, E = -70, Iapp = 0)
  formatphase!(Axis(F[14:18,3]), Axis(F[13,3]), args; res = res, tmax = 20000,
    ca0 = .55, x0 = .1, kwargs...)

##
  """  args = (xs = -4, cs = 0, gsyn = 0.01, E = -70, Iapp = 0)
    formatphase!(Axis(F[14:18,1]), Axis(F[13,1]), args; res = res, tmax = 20000)
    rowgap!(F.layout,3, Relative(.01))

    args = (xs = -4, cs = 0, gsyn = 0.0005, E = 30, Iapp = 0)
    formatphase!(Axis(F[14:18,2]), Axis(F[13,2]), args; res = res, tmax = 20000)

    args = (xs = 0, cs = -100, gsyn = 0.01, E = -70, Iapp = 0)
    formatphase!(Axis(F[14:18,3]), Axis(F[13,3]), args; res = res, tmax = 20000)
  """


  F
end

save("phase_examples_3row.jpeg", F)

"""f = Figure()
res = 30
args = (xs = -2.3, cs = 6, gsyn = 0, E = -70, Iapp = 0)
formatphase!(Axis(f[2:6,1]), Axis(f[1,1]), args; res = res, tmax = 900000, ca0 = 1.01, x0 = .75, cabeg = .5, caend = 1.5, xbeg = .4, xend = .9)
f

2

save("no_stable_basin_xs(2_4)_cs12.jpg", f)
"""
