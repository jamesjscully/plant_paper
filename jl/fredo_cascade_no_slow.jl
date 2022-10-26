using GLMakie, Colors, ColorSchemes
set_window_config!(framerate=3.0)
include("./plant.jl")
include("./compute_curves.jl")
include("./scan.jl")
include("./slow_transition_plot.jl")
begin
    α1 = .02
    β1 = .0011
    scale1 = α1/(α1+β1)
    α2 = .018
    β2 = .001
    scale2 = α2/(α2+β2)

    g12 = .011
    g21 = .03

    c1 = PlantCell(:c1, cs = -63, xs = -4, vout = [:s12], c0 = .4, v0 = 30, μc = 0, μx= 0, x0 = .72)
    c2 = PlantCell(:c2, cs = -10, xs = -4, vout = [:s21], c0 = 1.2, μc = 0, μx = 0, x0 = .5)
    s12 = LogisticSynapse(:s12, synout = [:c2], α = α1, β = β1, g = g12*scale1,
        ϵ = .00001, ♉ = -30, k = .1, E = 40, S0 = 0)
    s21 = LogisticSynapse(:s21, synout = [:c1], α = α2, β = β2, g = g21*scale2,
        ϵ = .00001, ♉ = -35, k = .1, E = -80, S0 = 0)

    net = Network([c1, c2, s12, s21])
    n = generate(net())

    prob = ODEProblem(n.f, n.u0, (0.,100000.))
    sol = solve(prob, BS3(), saveat = 2.)

    #smooth
    include("savitsky_golay.jl")
    sg1 = SavitzkyGolayFilter{250, 2}()
    sg2 = SavitzkyGolayFilter{100, 2}()
    #sol = smooth_trace(sol, net; f = sg1 ∘ sg1∘ sg1∘ sg2∘ sg2)
    v1 = [e.v_c1 for e in sol.u]
    v2 = [e.v_c2 for e in sol.u]
    s1 = [e.S_s12 for e in sol.u] #|> sgf |> sgf
    s2 = [e.S_s21 for e in sol.u]

    """:diverging_rainbow_bgymr_45_85_c67_n256v
    :rainbow1
    :jet1
    :gist_rainbow"""
    #push!(idx_bounds, 16000)
    F = Figure()
    t = sol.t
    ax = Axis(F[1,1])
    lines!(ax,t, v1)
    ax2 = Axis(F[2,1])
    lines!(ax2,t, v2)
    ax3 = Axis(F[3,1])
    lines!(ax3,t, s1)
    ax4 = Axis(F[4,1])
    lines!(ax4,t,s2)


    sol = smooth_trace(sol, net; f = sg1∘sg1)
    v1 = [e.v_c1 for e in sol.u]
    v2 = [e.v_c2 for e in sol.u]
    s1 = [e.S_s12 for e in sol.u] #|> sgf |> sgf
    s2 = [e.S_s21 for e in sol.u]
    h1 = [e.h_c1 for e in sol.u]
    h2 = [e.h_c2 for e in sol.u]
    idx_bounds = findsyncrits(net, s1  ,s2 )

    lines!(ax,t, v1)
    lines!(ax2,t, v2)
    lines!(ax3,t, s1)
    lines!(ax4,t,s2)
    is = 35:40
    scatter!(ax4, [t[i] for i in idx_bounds[is]], [s2[i] for i in idx_bounds[is]])

    F
end

F = network_cascade!(sol, idx_bounds[is];
    res = 15, cabeg = 0., caend = 1.3, xbeg = 0., N = 25)

save("fredo_traces_no_slow.jpeg", F)

lines(s1,s2,)
