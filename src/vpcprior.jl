l2 = 0.6931471805599453 # log(2)

function vpcprior(m::Float64, v::Float64, a::Float64, b::Float64, λ::Float64)
  d = digamma(a+b)-m*digamma(a)-(1-m)*digamma(b)

  z = log(λ)-
    λ*d^0.5+
    log(m*(1-m)*abs(trigamma(a)*m^2+trigamma(b)*(1-m)^2-trigamma(a+b)))-
    (l2+2*log(v)+0.5*log(d))

  return z
end

logdvdw(w::Float64) = w-2*(log(1+exp(w))+l2)

wpcprior(m::Float64, v::Float64, w::Float64, a::Float64, b::Float64, λ::Float64) = vpcprior(m, v, a, b, λ)+logdvdw(w)
