using Statistics
# this function is for getting dx and dc in the format required for streamplot.
# only use for streamplot.
function generate_fast_derivative_f(cs, xs, gsyn, E, Iapp)
    net = Network([FastCell(:c)])
    n = generate(net())
    u0 = n.uType([.2,.2,.2])
    fastf = n[1]
    function(c,x)
      prob = ODEProblem(fastf,u0,(0.,1500.),(x,c,gsyn,E, Iapp))
      sol = solve(prob, BS3(), saveat = 5.)
      cycles = 1 # minimum number of cycles to qualify as a spiker.
      cutoff = ceil(Int32, length(sol)/3)
      vtrace = sol[1,cutoff:end]
      idxs = [i for i in Iterators.drop(eachindex(vtrace),1) if |(
        vtrace[i]*vtrace[i-1] < 0 , vtrace[i] == 0)]
      dxavg = 0f0; dcavg = 0f0; avg = 0f0; freq = 0f0
      #if the number of zeros less than 2 then return the last value else return the trace between the last 2
      if length(idxs) < 2 * cycles + 1    #quiescent
          avg = last(vtrace)
          dxavg = dx(avg,x,xs)
          dcavg = dc(avg,x,c,cs)
      else  #spiker
          firstidx = length(idxs) % 2 == 0 ? 2 : 1
          tr = vtrace[idxs[firstidx]:idxs[end]]
          avg = convert(Float32, mean(tr))
          freq = convert(Float32, floor((length(idxs)-1)/2) /
            (sol.t[idxs[end]] - sol.t[idxs[firstidx]]))
          dxavg = mean(dx.(tr,x,xs))
          dcavg = mean(dc.(tr,x,c,cs))
      end
      return Point2(dcavg, dxavg)
    end
  end

  netm = Network([FastCell(:c)])
  nm = generate_ensemble(netm())
  function slowscan(Xsp, Csp, cs, xs, gsyn, E, Iapp)
      space = collect(Iterators.product(Xsp, Csp))
      res = length(Xsp)
      function output_func(sol, i)
          cycles = 2 # minimum number of cycles to qualify as a spiker.
          cutoff = ceil(Int32, length(sol)/5)
          vtrace = sol[1,cutoff:end]
          #get every zero by checking if there is a sign change between indices or if equal zero
          idxs = [i for i in Iterators.drop(eachindex(vtrace),1) if |(vtrace[i]*vtrace[i-1] < 0 , vtrace[i] == 0)]
          #calculate x and c from the index and their spaces
          x, c,=  space[i]
          dxavg = 0f0; dcavg = 0f0; avg = 0f0; freq = 0f0
          #if the number of zeros less than 2 then return the last value else return the trace between the last 2
          if length(idxs) < 2 * cycles + 1    #quiescent
              avg = last(vtrace)
              dxavg = dx(avg,x,xs)
              dcavg = dc(avg,x,c,cs)
          else  #spiker
              firstidx = length(idxs) % 2 == 0 ? 2 : 1
              tr = vtrace[idxs[firstidx]:idxs[end]]
              avg = convert(Float32, mean(tr))
              freq = convert(Float32, floor((length(idxs)-1)/2) / (sol.t[idxs[end]] - sol.t[idxs[firstidx]]))
              dxavg::Float32 = mean(dx.(tr,x,xs))
              dcavg::Float32 = mean(dc.(tr,x,c,cs))
          end
          return ((avg::Float32, freq::Float32, dxavg::Float32, dcavg::Float32),false)
      end

      function prob_funcm(prob,i,repeat)
          remake(prob,
              p = (space[i]..., gsyn, E, Iapp),
              u0 = rand(Float32, length(probm.u0))./10
          )
      end

      probm = ODEProblem{true, SciMLBase.FullSpecialize}(nm[1], nm.u0, (0f0,3000f0), (space[1], gsyn, E, Iapp))

      monteprob = EnsembleProblem(probm;
          prob_func = prob_funcm,
          output_func = output_func,
          reduction = (u, data, I) -> (push!(u, data), false),
          u_init = [],
          safetycopy = false)

      solscan = solve(monteprob,
          RK4(),
          EnsembleSerial(),
          batch_size = length(Xsp)*length(Csp),
          saveat = 15f0,
          trajectories=length(space))

          avg, freq, dxarr, dcarr = [collect(transpose(reshape(
            getindex.(solscan.u[1],i), res, res))) for i in 1:4]

  end
