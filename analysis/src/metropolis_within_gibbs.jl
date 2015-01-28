### This is not a generically coded Metropolis-within-Gibbs sampler, it is tailored to the problem at hand

function metropolis_within_gibbs(prior::ContinuousUnivariateDistribution,
                                 m::Vector{Float64},
                                 n::Matrix{Float64},
                                 x::Matrix{Float64},
                                 burnin::Int,
                                 nsteps::Int,
                                 Σ::Matrix{Float64}=eye(length(init)),
                                 thinning::Int=1)
  # Indices (i, j) refer to the i-th site (transcript) on the j-th cell
  # data_i = (mean_i, {n_ij: j}, {x_ij: j}), mean_i = x_i/n_i
  # Step 1: ∀ i draw a sample from v_i | ({p_ij: j}, data_i) using Metropolis sampling
  # Step 2: ∀ (i, j) draw a sample from p_ij | (v_i, data_i) which follows a Beta distribution

  nsites, ncells = size(n)
  @assert length(m) == nsites "length(m) must be equal to size(n, 1)"
  @assert size(x) == (nsites, ncells) "size(x) must be equal to size(n)"

  v = Array(Float64, nsites)
  p = Array(Float64, nsites, ncells)

  v[1] = rand(prior)
  a, b = beta_pars_from_mv(m, v[1])
  for j in 1:ncells
    p[1, j] = rand(Beta(a+x[1, j], b+n[1, j]-x[1, j]))
  end

  for i in 2:nsites
    a, b = beta_pars_from_mv(m, v[i-1])
    f(v::Float64) = logpdf(prior, v)+sum([logpdf(Beta(a+x[i, j], b+n[i, j]-x[i, j]), p[i-1, k]) for k in 1:ncells])
    v[i] = metropolis(f, v[i-1], burnin, nsteps, Σ, thinning)[end, :][1]

    for j in 1:ncells
      a, b = beta_pars_from_mv(m, v[i])
      p[i, j] = rand(Beta(a+x[i, j], b+n[i, j]-x[i, j]))
    end
  end

  return v, p
end
