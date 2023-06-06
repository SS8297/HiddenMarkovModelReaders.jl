####################################################################################################

"""

    setup(ɒ::M)
      where M <: Matrix{N}
      where N <: Number

# Description
Initiate hidden Markov model object.


See also: [`process!`](@ref), [`HMM`](@ref), [`HMMParams`](@ref)
"""
function setup(ɒ::M) where M <: Matrix{N} where N <: Number
  hmm = HMM([zeros(size(ɒ, 2), 1)], [zeros(size(ɒ, 1) + 1)], Array{Float64}(undef, 0))
  reset!(hmm)
  return hmm
end

####################################################################################################
"""

    setup(ɒ::M)
      where M <: Matrix{N}
      where N <: Number

# Description
Initiate hidden Markov model object.


See also: [`process!`](@ref), [`HMM`](@ref), [`HMMParams`](@ref)
"""
function setup(ɒ::A) where A ::Array{Float64, 3}
  hmm = HMM([zeros(ssize(ɒ, 2), size(ɒ, 3))], [zeros(size(ɒ, 1) + 1)], Array{Float64}(undef, 0))
  reset!(hmm)
  return hmm
end

####################################################################################################