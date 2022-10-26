using GLMakie
include("./plant.jl")
include("./compute_curves.jl")
include("./scan.jl")

function phase_plot!(ax, xs, cs, gsyn, E, Iapp; res = 100, stream=true, heat=true,
  traj = true, ca0 = .1, x0 = .1, xbeg = 0, xend = 1, cabeg = 0, caend = 2,
  tmax = 10000.)

  tmx = convert(Float64, tmax)

  limits!(ax, cabeg, caend, xbeg, xend)
  l = ax.limits[];

  #define search space and scan
  Xsp = range(l[2][1], stop = l[2][2], length = res)
  Csp = range(l[1][1], stop = l[1][2], length = res)
  avg, freq, dxarr, dcarr = slowscan(Xsp,Csp,cs,xs, gsyn, E, Iapp)

  # background heatmap of average v
  if heat
    heatmap!(ax, Csp, Xsp, avg, colormap = :coolwarm, colorrange = (-70, -20))
  end

  if stream
    f = generate_fast_derivative_f(cs, xs, gsyn, E, Iapp)
    streamplot!(ax, (c,x)-> Base.invokelatest(f,c,x),l[1][1] .. l[1][2],l[2][1] .. l[2][2],
      transparency = false, overdraw = true, arrow_size = 16, density = 1, linewidth = 1.5,
      maxsteps = 500, stepsize = .003, gridsize = (15, 30, 25),
      colormap = :Greys_9, colorrange = (0, 1))
  end

  #snic
  contour!(ax, Csp,Xsp,freq, levels = [0], color = :red, linewidth = 1)

  #Plot the trajectory
  if traj
    cell = PlantCell(:cell, xs = xs, cs = cs, gsyn = gsyn, E = E, Iapp = Iapp)
    net = Network([cell])
    n = generate(net())

    u0= n.uType([-60,0.,0.,x0,ca0])
    fullf = eval(n[2])
    prob = ODEProblem((u,p,t) -> Base.invokelatest(fullf,u,p,t),u0,(0.,tmx))
    sol = solve(prob, RK4())

    cidx = idx(net, :c_cell)
    xidx = idx(net, :x_cell)

    c2 = [e[cidx] for e in sol.u]
    x2 = [e[xidx] for e in sol.u]

    lines!(ax, c2, x2, color = :black, linewidth = 2)
  end

  #xavg nullcline using scan data
  contour!(ax, Csp,Xsp,dxarr, levels = [0], color = :darkblue,
  linewidth = 3)

  #xnullc analytic
  plot_x_nullcline_branches!(ax,xs,gsyn,E,Iapp; unstable_only = true, unstablecolor = :lightblue, linewidth = 3.2)

  # analytic ca nullcline
  #v = collect(-100:.1:30)
  #Ieq = get_fast_equilibrium(v, gsyn, E, Iapp)
  #c1, x1 = c_nullcline(v, cs, Ieq)
  #lines!(c1, x1, color = :green, linewidth = 2)

  # avg based ca nullcline
  contour!(ax, Csp,Xsp,dcarr, levels = [0], color = :grey, linewidth = 2)

  ax, sol

end

"""F = let F = Figure()
# xshift, cashift, gsyn, E, Iapp. See definition for optional kw arguments.
  ax1, sol1  = phase_plot!(Axis(F[2:6,1]),-4.,-100.,0.,-70,0; res = 100)
  hidespines!(ax1)
  hidedecorations!(ax1)

  ax2 = Axis(F[1,1])
  ax2.ylabel[] = "V(t)"
  lines!(ax2, sol1.t, sol1[1,:], color = :black)
  hidedecorations!(ax2, label = false)
  hidespines!(ax2)

  rowgap!(F.layout,1, Relative(.01))


  #ax2 = xnulc_fan_plot(Axis(F[1,2]))

  F
end"""

#save("phase_plot.jpeg", F)
##
