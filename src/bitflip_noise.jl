
function add_bitflip_noise(message, rate)
  n = length(message)
  noise = rand(Float64, n) .< rate
  noisy_message = message .⊻ noise
  return noisy_message
end
