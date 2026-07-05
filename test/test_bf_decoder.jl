@testitem "BitFlip Decoder" begin
using Test
using LDPCDecoders

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

  @test test_bf_decoder()

  # There is a low possibility of these tests failing
  @test test_bf_decoder_batch() < 0.005
end
