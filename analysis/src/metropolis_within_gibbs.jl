### This is not a generically coded Metropolis-within-Gibbs sampler, it is tailored to the problem at hand

function metropolis_within_gibbs(data::Dict{Symbol, Any},
                                 gibbs_prior::Dict{Symbol, Any},
                                 gibbs_init::Dict{Symbol, Any},
                                 gibbs_runner::Dict{Symbol, Any},
                                 metropolis_prior::Dict{Symbol, Any},
                                 metropolis_runner::Dict{Symbol, Any})
  nsites, ncells = size(data[:coverage])
  @assert length(data[:m]) == nsites "length(data[:m]) must be equal to size(data[:coverage], 1)"
  @assert size(data[:edited]) == (nsites, ncells) "size(data[:edited]) must be equal to size(data[:coverage])"


  chain_coords = (gibbs_runner[:burnin]+1):(gibbs_runner[:thinning]):(gibbs_runner[:nsteps])
  nsamples = length(chain_coords)

  previous_sample = Dict{Symbol, Any}(:v=>gibbs_init[:v], :p=>gibbs_init[:p])
  current_sample = Dict{Symbol, Any}(:v=>Array(Float64, nsites), :p=>Array(Float64, nsites, ncells))

  mcchain = Dict{Symbol, Any}(:v=>Array(Float64, nsamples, nsites), :p=>Array(Float64, nsamples, nsites, ncells))

  counter::Int = 1
  for k in 1:gibbs_runner[:nsteps]
    for i in 2:nsites
      a, b = beta_pars_from_mv(data[:m], previous_sample[:v][i])
      f(v::Float64) = logpdf(gibbs_prior[:v][i], v)+sum([logpdf(Beta(a+data[:edited][i, j], b+data[:coverage][i, j]-data[:edited][i, j]), previous_sample[:p][i, m]) for m in 1:ncells])
      current_sample[:v][i] = metropolis(f, previous_sample[:v][i], metropolis_runner[:burnin], metropolis_runner[:nsteps], metropolis_prior[:Σ], metropolis_runner[:thinning])[end, :][1]

      for j in 1:ncells
        a, b = beta_pars_from_mv(data[:m], current_sample[:v][i])
        current_sample[:p][i, j] = rand(Beta(a+data[:edited][i, j], b+data[:coverage][i, j]-data[:edited][i, j]))
      end
    end

    if in(k, chain_coords)
      mcchain[:v][counter, :] = current_sample[:v]
      mcchain[:p][counter, :, :] = current_sample[:p]
      counter += 1
    end

    copy!(previous_sample[:v], current_sample[:v])
    copy!(previous_sample[:p], current_sample[:p])
  end

  return mcchain
end
