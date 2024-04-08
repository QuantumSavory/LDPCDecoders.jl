using Test
using LDPCDecoders

@testset "test_bf_decoder.jl" begin

  """Test for bitflip decoder"""
  function test_bf_decoder()
    H = LDPCDecoders.parity_check_matrix(1000, 10, 9)
    per = 0.01
    err = rand(1000) .< per
    syn = (H * err) .% 2

    decoder = BitFlipDecoder(H, per, 100)
    guess, success = decode!(decoder, syn)

    return guess == err
  end

  """Test for batch bitflip decoder"""
  function test_bf_decoder_batch()
    H = LDPCDecoders.parity_check_matrix(1000, 10, 9)
    per = 0.01
    num_trials = 500
    errors = rand(1000,num_trials) .< per
    syndromes = H*errors .% 2
    decoder = BitFlipDecoder(H, per, 100)
    guesses, successes = batchdecode!(decoder, syndromes, zero(errors))
    actual_successes = [guesses[:,i]==errors[:,i] for i in 1:num_trials]
    ler = (num_trials - sum(actual_successes)) / num_trials
    return ler
  end

  """Test to verify old syndrome it decoder with new interface"""
  function test_deprecated_syndrome_it_decoder()
    H = LDPCDecoders.parity_check_matrix(1000, 10, 9)
    per = 0.01
    num_trials = 100
    errors = rand(Base.Float64, (1000,num_trials)) .< per
    syndromes = H*errors .% 2

    decoder = BitFlipDecoder(H, per, 100)
    guess_errors, success = batchdecode!(decoder, syndromes, zero(errors))
    actual_successes = 0
    matches_old_implementation = 0

    for i in axes(syndromes, 2)
      args = decoder.sparse_H, syndromes[:, i], decoder.max_iters, decoder.scratch.err, decoder.scratch.votes
      decoded_error, converged = LDPCDecoders.syndrome_it_decode(args...)
      actual_successes += decoded_error == errors[:, i]
      
      # Small change in new implementation, most of the guesses should match
      matches_old_implementation += decoded_error == guess_errors[:, i]
    end

    @assert matches_old_implementation >= 0.90 * num_trials
    ler = (num_trials - sum(success)) / num_trials
    return ler
  end

  @test test_bf_decoder()

  # There is a low possibility of these tests failing
  @test test_bf_decoder_batch() < 0.005
  @test test_deprecated_syndrome_it_decoder() < 0.005
end
