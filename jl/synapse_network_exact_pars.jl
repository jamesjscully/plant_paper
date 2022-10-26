using OrdinaryDiffEq, GLMakie
using Statistics:mean

include("./plant.jl")

α1 = .02
β1 = .0011
scale1 = α1/(α1+β1)
α2 = .018
β2 = .001
scale2 = α2/(α2+β2)

g12 = .011
g21 = .03

## get coupling function when c2 is postsynaptic cell s12
begin
    varr = Array{Float64}(undef,0)
    const sp = 0. : .05 : 1
    for s in sp
        gsyn = s*g12*scale1
        c = PlantCell(:c, cs = -10, xs = -4, vout = [:s21],
            c0 = 1.2, μc = 0, μx = 0, x0 = .5, E = 40, gsyn = gsyn)
        net = Network([c])
        n = generate(net())
        prob = ODEProblem(n.f, n.u0, (0.,10000.))
        sol = solve(prob, BS3(), saveat = 1.0)
        v = [e.v_c for e in sol.u]
        push!(varr, mean(v[ceil(Int64,length(v)/2):end]))
    end

    const _varr = deepcopy(varr)
    function vpre1(s)
        N = length(_varr)
        i = s==0 ? 1 : ceil(Int64, s*(N-1))
        x0 = sp[i]
        x1 = sp[i+1]
        y0 = _varr[i]
        y1 = _varr[i+1]
        return (y0*(x1-s) + y1*(s-x0))/(x1-x0)
    end


    ## get coupling function when c1 is postsynaptic cell s21
    varr2 = Array{Float64}(undef,0)

    for s in sp
        gsyn = s*g21*scale2
        c = PlantCell(:c, cs = -63, xs = -4, vout = [:s],
            c0 = .4, v0 = 30, μc = 0, μx= 0, x0 = .72, E = -80, gsyn = gsyn)
        net = Network([c])
        n = generate(net())
        prob = ODEProblem(n.f, n.u0, (0.,10000.))
        sol = solve(prob, BS3(), saveat = 1.0)
        v = [e.v_c for e in sol.u]
        push!(varr2, mean(v[ceil(Int64,length(v)/2):end]))
    end

    const _varr2 = deepcopy(varr2)
    function vpre2(s)
        N = length(_varr2)
        i = s==0 ? 1 : ceil(Int64, s*(N-1))
        x0 = sp[i]
        x1 = sp[i+1]
        y0 = _varr2[i]
        y1 = _varr2[i+1]
        return (y0*(x1-s) + y1*(s-x0))/(x1-x0)
    end


    u0 = [0.2,0.8]
    p = (
        α1 = α1,
        β1 = β1,
        ϵ1 = .00001,
        ♉1 = -30.,
        k1 = .1,
        α2 = α2,
        β2 = β2,
        ϵ2 = .00001,
        ♉2 = -35.,
        k2 = .1,
    )

    tspan = (0., 100000.)

    function synaptic_net(du,u,p,t)
        α1,β1,ϵ1,♉1,k1,α2,β2,ϵ2,♉2,k2 = p
        S1,S2 = u
        du[1] = α1*S1*(1-S1)/(1+exp(♉1-vpre1(S2))/k1) - β1*(S1 -ϵ1)
        du[2] = α2*S2*(1-S2)/(1+exp(♉2-vpre2(S1))/k2) - β2*(S2 -ϵ2)
    end

    prob = ODEProblem(synaptic_net,u0,tspan,p)

    sol = solve(prob, RK4())
    lines(sol[1,:], sol[2,:])
end
