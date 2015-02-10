logit(v::Float64) = log(v/(1-v))

inv_logit(w::Float64) = exp(w)/(1+exp(w))
