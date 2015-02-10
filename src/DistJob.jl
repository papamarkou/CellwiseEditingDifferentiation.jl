### The DistJob type is defined as a hint for a future fully-treated generic implementation of the Gibbs algorithm
### In fact the existing job types in Lora can accommodate DistJob after minor modifications of DistJob

type DistJob
  target::Distribution
  #nsamples::Int
end

Job = Union(MCJob, DistJob)

run(j::DistJob) = rand(j.target)
