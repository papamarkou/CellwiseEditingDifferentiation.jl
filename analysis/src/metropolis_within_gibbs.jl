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

  mcchain = Dict{Symbol, Any}()
  mcchain[:v] = Array(Float64, nsites, nsamples)
  mcchain[:p] = Array(Float64, nsites, ncells, nsamples)

  counter::Int = 1
  for i in 2:nsites
    a, b = beta_pars_from_mv(data[:m][i], mcchain[:v][i-1])
    f(v::Float64) = logpdf(gibbs_prior[:v][i], v)+sum([logpdf(Beta(a+data[:edited][i, j], b+data[:coverage][i, j]-data[:edited][i, j]), p[i-1, k]) for k in 1:ncells])
    v[i] = metropolis(f, v[i-1], metropolis_runner[:burnin], metropolis_runner[:nsteps], metropolis_prior[:Σ], metropolis_runner[:thinning])[end, :][1]

    for j in 1:ncells
      a, b = beta_pars_from_mv(data[:m], v[i])
      p[i, j] = rand(Beta(a+data[:edited][i, j], b+data[:coverage][i, j]-data[:edited][i, j]))
    end
  end

  return mcchain
end
