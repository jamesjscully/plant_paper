# this file defines the models for the full plant model and fast subsystem.
# the package CoupledODETools is not documented, but you should not need to make
# changes.

using CoupledODETools
using OrdinaryDiffEq

# helper functions
r(v)= 127v/105 +8265/105
function m∞(v)
    amm = @. 0.1*(50-r(v))/(exp((50-r(v))/10)-1)
    bmm = @. 4*exp((25-r(v))/18)
    m∞ = @. amm/(amm+bmm)
end
H(x) = x>0 ? 1 : 0

@Component PlantCell begin
    Dv = gna*^(m∞(v),3.0)*h*(30-v)+gn*^(n,4.0)*(-75 -v)+.01x*(30-v)+.03c/(.5+c)*(-75-v)+gleak*(-40-v) + gsyn*(E-v) + Iapp + Iext , 0
    Dh = μh*((1-h)*(.07*exp((25-r(v))/20))-h*(1/(1+exp((55-r(v))/10)))),0
    Dn = μn*((1-n)*(.01*(55-r(v)))/(exp((55-r(v))/10)-1)-n*0.125*exp((45 - r(v))/80)),.0
    Dx = μx*((1/(exp(.15*(-50 + xs -v))+1))-x), .5
    Dc = μc*(.0085*x*(140 - v +cs)-c), .2
    vout = vpre -> ~v
    voutf = vfasc -> ~v
end begin
    μh = .08
    μn = .08
    gna = 4
    gn = .3
    xs = -4
    cs = -60
    gleak = .003
    μx = .012
    μc = .00025
    gsyn = 0
    E = 0
    Iapp = 0
    Iext = 0
end

@Component FastCell begin
    Dv = gna*^(m∞(v),3.0)*h*(30-v)+gn*^(n,4.0)*(-75 -v)+.01p[1]*(30-v)+
        .03p[2]/(.5+p[2])*(-75-v)+gleak*(-40-v) + p[3]*(p[4]- v) + p[5],0
    Dh = μh*((1-h)*(.07*exp((25-r(v))/20))-h*(1/(1+exp((55-r(v))/10)))),0
    Dn = μn*((1-n)*(.01*(55-r(v)))/(exp((55-r(v))/10)-1)-
        n*0.125*exp((45 - r(v))/80)),.0
end begin
    μh = .08
    μn = .08
    gna = 4
    gn = .3
    gleak = .003
end

@Component Diode begin
    synout = Iext -> g*(vpre-v)
end begin
    g = .1
end

@Component DynamicSynapse begin
    DS = α*(1-S)/(1+exp(♉s-vpre)/k) - β*S, 0
    DM = 1/(1+exp(♉m-vpre) - M)/τ, 0
    synout = Iext -> g*~S*~M*(E-v)
    synout2 = Iext -> g*~S^2*~M*(E-v)
end begin
    α = .1
    β = .1
    ♉s = 0
    ♉m = -20
    k = 20
    E = -60
    g = 0.1
    τ = 1000
end
@Component AlphaSynapse begin
    DS = α*(1-S)/(1+exp(♉-vpre)/k) - β*S, 0
    synout = Iext -> g*~S*(E-v)
end begin
    α = .1
    β = .1
    ♉ = 0
    k = 20
    E = -60
    g = 0.1
end

@Component LogisticSynapse begin
    DS = α*S*(1-S)/(1+exp(♉-vpre)/k) - β*(S -ϵ), .002
    synout = Iext -> g*~S*(E-v)
end begin
        α = .1
        β = .1
        ϵ = .001
        ♉ = -15
        k = 20
        E = -60
        g = .1
end

# some utilities for calculating burst period and such.
function bperiod(sol,i)
    #get every zero by checking if there is a sign change between indices or if equal zero
    vtrace = sol[idx,:]
    idxs = filter(2:length(vtrace)) do i
        (vtrace[i-1] > 0) & (vtrace[i] < 0)
    end
    times = diff([sol.t[i] for i in idxs]) |> sort
    if length(times) < 3 return 0.0 end
    biggest_increase = mapreduce( (x,y) -> x.increase > y.increase ? x : y, collect(2:length(times))) do i
        (increase = (times[i]-times[i-1])/times[i], idx = i)
    end
    mean(times[biggest_increase.idx:end])/500
end
function hcomeasure(data::Array{T}, time::Array{T}, thresh = 0) where T
    last_up2down = time[1]
    dtimes = T[]
    #calculate the time between roots when f`<0
    for i = 1:length(data)-1
        if ((data[i+1]-thresh)*(data[i]-thresh) < -.001) & (data[i+1] < data[i])
            push!(dtimes,time[i+1]-last_up2down)
            last_up2down = time[i+1]
        end
    end
    #ratio of minimum 3 to maximum
    if length(dtimes) < 11 return 0.0
    else
        sorted = sort(dtimes) # dump outliers. Only one for big
        mode = ceil(Int32, length(sorted)/2)
        midn = sorted[mode-5:mode+5]; lastn = sorted[end-3:end]
        return ratio = mean(lastn)/mean(midn)
    end
end

# used to get the index of a certain variable.
idxs(net) = keys(generate(net()).u0)
idx(net, sym) = findfirst(x -> x == sym, idxs(net))

@inline function dx(v,x,sft)
    .012*((1/(exp(0.15*(-v-50+sft)) + 1))-x)
end

@inline function dc(v,x,c,sft)
    0.00025*(0.0085*x*(140-v+sft)-c)
end
