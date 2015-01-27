### Metropolis sampler

# q: proposal density
# init: initial value of parameters
function metropolis(q, init, burnin, thinning, nsteps)
  r = (burnin+1):thinning:nsteps
  npars = length(init)
  nsamples = length(r)

  mcchain = Matrix(Float64, nsamples, npars)
  pars = Vector(Float64, npars)

  for i in r
    pars = rand(q)
  end
end
