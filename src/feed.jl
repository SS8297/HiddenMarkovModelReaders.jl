####################################################################################################

"Feed on hidden Markov model object."
function feed!(∫::HMM, frame::I, ɒ::A{F, 2}; params::HMMParams) where I <: Int64 where A <: Array where F <: Float64
  for ι ∈ eachindex(∫.model)
    plus = euclDist(∫.data[ι], ɒ[frame, :])
    for ο ∈ eachindex(∫.model)
      lpen = copy(params.penalty)
      if (ο == ι) lpen = 0 end
      if ∫.model[ο][frame + 1] < 0. || ∫.model[ι][frame] + plus + lpen < ∫.model[ο][frame + 1]
        ∫.model[ο][frame + 1] = ∫.model[ι][frame] + plus + lpen
      end
    end
  end
end

####################################################################################################
