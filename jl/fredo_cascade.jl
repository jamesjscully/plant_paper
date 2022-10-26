using GLMakie, Colors
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

    c1 = PlantCell(:c1, cs = -63, xs = -4, vout = [:s12], c0 = .5, v0 = 30)
    c2 = PlantCell(:c2, cs = -10, xs = -4, vout = [:s21], c0 = 1.2)
    s12 = LogisticSynapse(:s12, synout = [:c2], α = α1, β = β1, g = g12*scale1,
        ϵ = .00001, ♉ = -30, k = .1, E = 40, S0 = 0)
    s21 = LogisticSynapse(:s21, synout = [:c1], α = α2, β = β2, g = g21*scale2,
        ϵ = .00001, ♉ = -35, k = .1, E = -80, S0 = 0)

    net = Network([c1, c2, s12, s21])
    n = generate(net())

    prob = ODEProblem(n.f, n.u0, (0.,100000.))
    sol = solve(prob, BS3())
end

#smooth

v1 = [e.v_c1 for e in sol.u]
v2 = [e.v_c2 for e in sol.u]
s1 = [e.S_s12 for e in sol.u] #|> sgf |> sgf
s2 = [e.S_s21 for e in sol.u] #|> sgf |> sgf
# plot the traces with index on bottom
let
    F = Figure()
    ax = Axis(F[1,1])
    lines!(ax, v1)
    ax2 = Axis(F[2,1])
    lines!(ax2, v2)
    ax3 = Axis(F[3,1])
    lines!(ax3, s1)
    ax4 = Axis(F[4,1])
    lines!(ax4, s2)
    F
end

# partition the trace and smooth it

let
    beg = 3850
    fin = 4300
    F = Figure()
    ax = Axis(F[1,1])
    lines!(ax, sol.t[beg:fin], v1[beg:fin])
    ax2 = Axis(F[2,1])
    lines!(ax2, sol.t[beg:fin], v2[beg:fin])
    ax3 = Axis(F[3,1])
    lines!(ax3, sol.t[beg:fin], s1[beg:fin])
    ax4 = Axis(F[4,1])
    lines!(ax4, sol.t[beg:fin], s2[beg:fin])
    F
end
let
    beg = 4300
    fin = 5020
    F = Figure()
    ax = Axis(F[1,1])
    lines!(ax, sol.t[beg:fin], v1[beg:fin])
    ax2 = Axis(F[2,1])
    lines!(ax2, sol.t[beg:fin], v2[beg:fin])
    ax3 = Axis(F[3,1])
    lines!(ax3, sol.t[beg:fin], s1[beg:fin])
    ax4 = Axis(F[4,1])
    lines!(ax4, sol.t[beg:fin], s2[beg:fin])
    F
end
let
    beg = 5020
    fin = 5200
    F = Figure()
    ax = Axis(F[1,1])
    lines!(ax, sol.t[beg:fin], v1[beg:fin])
    ax2 = Axis(F[2,1])
    lines!(ax2, sol.t[beg:fin], v2[beg:fin])
    ax3 = Axis(F[3,1])
    lines!(ax3, sol.t[beg:fin], s1[beg:fin])
    ax4 = Axis(F[4,1])
    lines!(ax4, sol.t[beg:fin], s2[beg:fin])
    F
end
let
    beg = 5200
    fin = 5250
    F = Figure()
    ax = Axis(F[1,1])
    lines!(ax, sol.t[beg:fin], v1[beg:fin])
    ax2 = Axis(F[2,1])
    lines!(ax2, sol.t[beg:fin], v2[beg:fin])
    ax3 = Axis(F[3,1])
    lines!(ax3, sol.t[beg:fin], s1[beg:fin])
    ax4 = Axis(F[4,1])
    lines!(ax4, sol.t[beg:fin], s2[beg:fin])
    F
end

idx_bounds = (3850,5020,5250)
#idx_bounds = (3850,4300,5020,5200,5250)


F = network_cascade!(Figure(resolution = (1800,1800)), sol, idx_bounds;
    res = 100, cabeg = .2, caend = 1.5, xbeg = .5, N = 6)


save("fredo_cascade.jpeg", F)
