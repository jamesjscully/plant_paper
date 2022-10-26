# this file defines the computations for snic, xdot=0, etc.
"""
this function returns the internal current of a steady state given v with
  optional kw args (after ;): gsyn = 0, E = -70, Iapp = 0)
"""
function get_fast_equilibrium(v, gsyn, E, Iapp)
  r = @. 127v/105 + 8265/105

  gna = 4
  amm = @. 0.1*(50-r)/(exp((50-r)/10)-1)
  bmm = @. 4*exp((25-r)/18)
  m∞ = @. amm/(amm+bmm) #@. 1/(exp((-22-v)/8)+1)
  _h1 = @. .07*exp((25-r)/20)
  _h2 = @. 1/(1+exp((55-r)/10))
  h = @. _h1 / (_h1 +_h2)
  INa = @. gna*^(m∞,3.0)*h*(30-v)

  gn = .3
  _n1 = @. .01*(55-r)/(exp((55-r)/10)-1)
  _n2 = @. 0.125*exp((45 - r)/80)
  n= @. _n1 / (_n1 +_n2)
  In = @. gn*^(n,4.0)*(-75 -v)

  gleak = .003
  Il = @. gleak*(-40-v)

  Isyn = @. gsyn*(E-v)

  return @. INa + In + Il + Isyn + Iapp
end

#analytic-ish curves of equilibrium
function x_nullcline(v, xs, Ieq)
  x = @. 1/(exp(.15*(-50 + xs -v))+1)
  Ix = @. .01x*(30-v)
  q = Ieq + Ix
  c = @. (-.015v-1.125)/(-.03v+q-2.25)-1/2
  c = [ e>0. ? e : 0. for e in c]
  return (c,x)
end

function c_nullcline(v, cs, Ieq)
  l = @. (-75-v)
  m = @. .01*(30-v)
  d = @. .0085*(140-v+cs)
  p = @. m/d

  a = @. -.5*Ieq
  b = @. .5p+.03l-Ieq
  c = @. Complex(p)

  c1 = @. (-b+sqrt(b^2.0 -4*a*c))/(2a)
  c2 = @. (-b-sqrt(b^2.0 -4*a*c))/(2a)
  x1 = @. c1/d
  x2 = @. c2/d
  # first branch is not used
  """
  c3 = Float64[]; x3 = Float64[]
  for i in eachindex(c1)
    c = real(c1[i])
    x = real(x1[i])
    if (isreal(c1) && isreal(x1) && (c>0) && (x<1) && (x>0))
      push!(c3,c); push!(x3,x)
    end
  end"""

  c4 = Float64[]; x4 = Float64[]
  for i in eachindex(c2)
    c = real(c2[i])
    x = real(x2[i])
    if (isreal(c2) && isreal(x2) && (c>0) && (x<1) && (x>0))
      push!(c4,c); push!(x4,x)
    end
  end
  return (c4, x4) # (c3, x3, c4, x4)
end

function x_nullcline_branches(v,xs,Ieq; upper = true, unstable_only = false)
  c,x = x_nullcline(v,xs,Ieq)
  #dx/dc
  slopes = diff(x)./diff(c)
  # find indexes when dx/dc changes sign
  idxs = [i for i in Iterators.drop(eachindex(slopes),1) if
    (slopes[i]*slopes[i-1] < 0) || (slopes[i] == 0)]
    if length(idxs) > 2
      fin =  c[idxs[3]] > 5 ? 2 : min(4,length(idxs)-1)
    else
      fin = length(idxs)-1
    end

  # get rid of irrelevant unstable upper branch
  if !upper
    if length(idxs) > 2 fin = min(3,fin) end
  end

  x_b = [x[idxs[j]:idxs[j+1]] for j = 1 : fin]
  c_b = [c[idxs[j]:idxs[j+1]] for j = 1 : fin]

  if unstable_only
    x_b = [x[idxs[j]:idxs[j+1]] for j = 2 : 2 : fin]
    c_b = [c[idxs[j]:idxs[j+1]] for j = 2 : 2 : fin]
  end

  return (c_b, x_b)
end

function plot_x_nullcline_branches!(ax,xs, gsyn, E, Iapp;
   stablecolor = :darkblue, unstablecolor = :darkblue, upper = true,
   unstable_only = false, kwargs... )
   v = -100:.1:30
   Ieq = get_fast_equilibrium(v, gsyn, E, Iapp)
   cl,xl = x_nullcline_branches(v,xs,Ieq; upper = upper,
      unstable_only = unstable_only)
   for i in eachindex(cl)
     c = cl[i]; x = xl[i]
     j= unstable_only ? 1 : 0
     style = (i+j%2)==1 ? :solid : :dot
     color = (i+j%2)==1 ? stablecolor : unstablecolor
     lines!(ax,c,x; linestyle = style, color = color, kwargs...)
   end
end
