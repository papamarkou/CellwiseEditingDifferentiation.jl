### Metropolis sampler

# f: log of target
# init: initial value of parameters
# burnin: burn-in length
# nsteps: total number of Monte Carlo iterations including burn-in length
# Σ: covariance matrix of normal proposal density
# thinning: thinning of returned Monte Carlo chain
function metropolis(f::Function,
                    init::Vector{Float64},
                    burnin::Int,
                    nsteps::Int,
                    Σ::Matrix{Float64}=eye(length(init)),
                    thinning::Int=1)
  npars = length(init)

  chain_coords = (burnin+1):thinning:nsteps
  nsamples = length(chain_coords)

  current_sample = copy(init)

  mcchain = Matrix(Float64, nsamples, npars)

  i::Int = 1
  for j in 1:nsteps
    proposed_sample = rand(MvNormal(current_sample, Σ))

    ratio = f(proposed_sample)-f(current_sample)

    if ratio > 0 || (ratio > log(rand()))
      current_sample = copy(proposed_sample)
    end

    if in(j, chain_coords)
      mcchain[i, :] = current_sample
      i += 1
    end
  end
end
