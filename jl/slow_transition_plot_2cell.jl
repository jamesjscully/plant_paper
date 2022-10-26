

# generate a trajectory

getpar(component, sym) = [x.val for x in component.pars if x.sym == sym][1]

function smooth_trace(sol,net; f = x -> x, vthresh = -35, tthresh = 1000 )
    cells = filter(x -> x isa PlantCell, net.components)
    segs = []
    u = hcat(Array.(sol.u)...)
    for c in cells
        v_name = Symbol("v_", c.name)
        v = getindex.(sol.u, v_name)
        upixs = filter(collect(1:length(v)-1)) do i
            ((v[i]-vthresh)*(v[i+1]-vthresh) < 0.) && (v[i] < v[i+1])
        end

        cellvars = Symbol.(filter(endswith(string(c.name)),
            string.([keys(sol.u[1])...])))
        synvars = [Symbol("S_",coup.targets[1]) for coup in c.couplings if
            coup.flag == :vpre]
        vars = vcat(cellvars,synvars)
        varixs = [idx(net,var) for var in vars]
        lastdt = 0
        for ix in 1:length(upixs)-1
            i = upixs[ix]; j = upixs[ix+1]
            dt = sol.t[j]-sol.t[i]
            ixrange = if dt < tthresh
                i:j
            else
                fini = findfirst(x -> x > (sol.t[i]+lastdt) , sol.t)
                if isnothing(fini)
                    i:length(sol.t)
                else
                    i:fini
                end
            end
            for vari in varixs
                u[vari, ixrange] = fill(mean(view(sol,vari,ixrange)), length(ixrange))
            end
            lastdt = dt
        end
    end
    for i = 1:length(sol.u[1])
        u[i,:] = f(view(u,i,:))
    end
    (u = typeof(sol.u[1]).(collect(eachcol(u))), t = sol.t)
end

function findsynmaxima(net, s...)
    tixs = sort(vcat([findall(i ->
        (diff(x)[i+1] < 0)&&(diff(x)[i] > 0 )&&(x[i]> mean(x)),
        1:length(x)-2) for x in s]...))
end
function findsyncrits(net, s...)
    tixs = sort(vcat([findall(i ->
        (diff(x)[i+1]*diff(x)[i]<0),
        1:length(x)-2) for x in s]...))
end



function network_cascade!(sol, spikesol, idx_bounds;
    res = 100, xbeg = 0, xend = 1, cabeg = 0, caend = 2.5, N=10)

    fig = Figure(resolution = (1400*(length(idx_bounds)-1),3000))
    set_theme!(fontsize = 35)

    cells = filter(x -> x isa PlantCell, net.components)
    v = [[getproperty(e, Symbol("v_",cell.name)) for e in sol.u] for cell in cells]
    c = [[getproperty(e, Symbol("c_",cell.name)) for e in sol.u] for cell in cells]
    x = [[getproperty(e, Symbol("x_",cell.name)) for e in sol.u] for cell in cells]
    spikev = [[getproperty(e, Symbol("v_",cell.name)) for e in spikesol.u] for cell in cells]
    synapses = vcat([filter(x -> (x isa LogisticSynapse)&(cell.name in x.couplings[1].targets),
            net.components)
        for cell in cells]...)
    g = [getpar(syn,:g) for syn in synapses]
    s = [[getproperty(e, Symbol("S_",syn.name)) for e in sol.u] for syn in synapses]
    spikes = [[getproperty(e, Symbol("S_",syn.name)) for e in spikesol.u] for syn in synapses]

    """c1 = colorant"red"
    c2 = colorant"blue"
    c3 = colorant"lightsalmon"
    c4 = colorant"lightblue"
    colors = range(c1, stop = c2, length = N)
    colors2 = range(c3, stop = c4, length = N)"""

    colors = get(colorschemes[:coolwarm], collect(1:N), :extrema)
    colors2 = colors # get(colorschemes[:diverging_rainbow_bgymr_45_85_c67_n256
    #], collect(1:N), :extrema)

    tt = (sol.t .- sol.t[idx_bounds[1]])./1000
    for i=1:length(idx_bounds)-1

        tpoints = range(sol.t[idx_bounds[i]],
            stop = sol.t[idx_bounds[i+1]], length = N+1)

        ixs = [findfirst(>(t), sol.t) for t in tpoints]

        tmidpoints = diff(collect(tpoints)) .+ tpoints[1:end-1]
        ixsmid = [findfirst(>(t), sol.t) for t in tmidpoints]
        gsynarr = [getindex(s[i].*g[i],ixsmid) for i in eachindex(s)]


        ax1 = Axis(fig[7,i],
            xticks = WilkinsonTicks(2),yticks = WilkinsonTicks(2),
            backgroundcolor = :grey)
        ax1.ylabel[] = L"S_{pre}(t)"
        ax1.xlabel[] = L"t\: [s]"
        ylims!(ax1,0,maximum(maximum.(getindex.(
            s, Ref(idx_bounds[1]: idx_bounds[end])))))
        xlims!(ax1, tt[ixs[1]], tt[ixs[end]])

        lines!(ax1,tt[i:i+1],s[1][i:i+1])
        lines!(ax1,tt[i:i+1],s[2][i:i+1])

        ax2 = Axis(fig[6,i], yticks = WilkinsonTicks(2),
            xticklabelsvisible = false, xticksvisible = false,
            backgroundcolor = :grey)
        ylims!(ax2,-80,30)
        xlims!(ax2, sol.t[ixs[1]]./1000, sol.t[ixs[end]]./1000)
        ax2.ylabel[] = L"V(t)\: [mV]"

        lines!(ax1,tt[i:i+1],v[1][i:i+1])
        lines!(ax1,tt[i:i+1],v[2][i:i+1])

        vm = minimum(vcat(v[1],v[2]))
        col = [:black, :white]
        for (j,cell) in enumerate(cells)

            ax3 = Axis(fig[1+(j-1)*7:5+(j-1)*7,i],
                xticks = WilkinsonTicks(2),yticks = WilkinsonTicks(2), backgroundcolor = :grey)
            ax3.ylabel[] = "x variable"
            ax3.xlabel[] = "[Ca]"
            limits!(ax3, cabeg, caend, xbeg, xend)
            l = ax3.limits[];

            Xsp = range(l[2][1], stop = l[2][2], length = res)
            Csp = range(l[1][1], stop = l[1][2], length = res)

            xs = getpar(cell, :xs)
            cs = getpar(cell, :cs)
            E = getpar(synapses[j], :E)
            Iapp = getpar(cell, :Iapp)

            """for ax in [ax1,ax2,ax3]
                ax.backgroundcolor[] = :grey
            end"""
            lines!(ax3, c[j][idx_bounds[1]:13071],
                x[j][idx_bounds[1]:13071], color = col[j], linewidth = 2)
            #println(maximum(gsynarr[j]))
            for (k, gsyn) in reverse(collect(enumerate(gsynarr[j])))
                avg, freq, dxarr, dcarr = slowscan(Xsp,Csp,cs,xs, gsyn, E, Iapp)
                #contour!(ax3, Csp, Xsp, freq, levels = [0f0],
                #    color = colors2[k], linewidth = .5)
                contour!(ax3, Csp, Xsp, dxarr, levels = [0f0],
                    color = colors[k], linewidth = 1)
                contour!(ax3, Csp, Xsp, dcarr, levels = [0f0],
                    color = :darkgrey, linewidth = 1)
                plot_x_nullcline_branches!(ax3, xs, gsyn, E, Iapp;
                    stablecolor = colors2[k], unstablecolor = colors2[k], upper = false, linewidth = 6)
            end

            for k in eachindex(gsynarr[j])
                lines!(ax3, c[j][ixs[k]:ixs[k+1]], x[j][ixs[k]:ixs[k+1]],
                    color = col[j], linewidth = 23)
                lines!(ax3, c[j][ixs[k]:ixs[k+1]], x[j][ixs[k]:ixs[k+1]],
                    color = colors[k], linewidth = 15)
                lines!(ax1, tt[ixs[k]:ixs[k+1]], spikes[j][ixs[k]:ixs[k+1]],
                    color = col[j], linewidth = .5)
                lines!(ax2, tt[ixs[k]:ixs[k+1]], spikev[j][ixs[k]:ixs[k+1]],
                    color = col[j], linewidth = .5)
            end
        end
        for j in eachindex(cells)
            for k in eachindex(gsynarr[1])
                lines!(ax1, tt[ixs[k]:ixs[k+1]], s[j][ixs[k]:ixs[k+1]],
                    color = col[j], linewidth = 7)
                lines!(ax2, tt[ixs[k]:ixs[k+1]], v[j][ixs[k]:ixs[k+1]],
                    color = col[j], linewidth = 7)
            end
        end
        for k in eachindex(gsynarr[1])
            t = tt[ixs[k]:ixs[k+1]]
            lines!(ax1, t, zeros(length(t)).- .1,
                color = colors[k], linewidth = 30)
            """lines!(ax2, t, fill(vm,length(t)),
                color = colors[k], linewidth = 15)"""
        end
        ylims!(ax1,(-.1,.5))
    end
    fig
end
