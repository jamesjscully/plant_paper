using GLMakie, ColorSchemes
set_window_config!(framerate=3.0)
include("./plant.jl")
include("./compute_curves.jl")
include("./scan.jl")
include("./slow_transition_plot_2_cell.jl")

begin
    cpars = (cs = -50,xs = -4)
    spars = (α = .05, β = .0051, g = .033, ϵ = .001, ♉ = -20, k = .1, E = -80, S0 = 0)

    c1 = PlantCell(:c1; vout = [:s12], c0 = .65, cpars...)
    c2 = PlantCell(:c2; vout = [:s21], c0 = .7, cpars...)

    s12 = LogisticSynapse(:s12; synout = [:c2], spars...)
    s21 = LogisticSynapse(:s21; synout = [:c1], spars...)

    net = Network([c1, c2, s12, s21])
    n = generate(net())

    prob = ODEProblem(n.f, n.u0, (0.,100000.))
    sol = solve(prob, BS3(), adaptive = false, dt = 2)

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


    sol2 = smooth_trace(sol, net; f = sg1∘sg1)
    v1 = [e.v_c1 for e in sol2.u]
    v2 = [e.v_c2 for e in sol2.u]
    s1 = [e.S_s12 for e in sol2.u] #|> sgf |> sgf
    s2 = [e.S_s21 for e in sol2.u]
    idx_bounds = findsynmaxima(net,s1 |> sg1 |> sg1 ,s2 |> sg1|> sg1)

    lines!(ax,t, v1)
    lines!(ax2,t, v2)
    lines!(ax3,t, s1)
    lines!(ax4,t,s2)
    scatter!(ax4, [t[i] for i in idx_bounds], [s2[i] for i in idx_bounds])

    F
end


F = network_cascade!(sol2,sol, idx_bounds[2:3];
    res = 15, cabeg = 0., caend = 1.3, xbeg = 0., N = 25)
set_theme!()

save("hco_cascade_tonic_stable_2.jpeg", F)
