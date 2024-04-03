using Test
using LDPCDecoders

@testset "test_bp_decoder.jl" begin

  """Test for belief propagation decoder"""
  function test_bp_decoder()
    H = LDPCDecoders.parity_check_matrix(1000, 10, 9)
    per = 0.01
    err = rand(1000) .< per
    syn = (H * err) .% 2

    bpd = BeliefPropagationDecoder(H, per, 100)
    guess, success = decode!(bpd, syn)

    return guess == err
  end

  """Test for batch belief propagation decoder"""
  function test_bp_decoder_batch()
    H = LDPCDecoders.parity_check_matrix(1000, 10, 9)
    per = 0.01
    num_trials = 100
    errors = rand(1000,num_trials) .< per
    syndromes = H*errors .% 2
    bpd = BeliefPropagationDecoder(H, per, 100)
    guesses, successes = batchdecode!(bpd, syndromes, zero(errors))
    actual_successes = [guesses[:,i]==errors[:,i] for i in 1:num_trials]
    ler = (num_trials - sum(actual_successes)) / num_trials
    return ler
  end

  """Test for batch belief propagation decoder"""
  function test_deprecated_syndrome_decoder()
    H = LDPCDecoders.parity_check_matrix(1000, 10, 9)
    per = 0.01
    num_trials = 100
    errors = rand(Base.Float64, (1000,num_trials)) .< per
    syndromes = H*errors .% 2

    decoder = BeliefPropagationDecoder(H, per, 100)
    guess_errors, success = batchdecode!(decoder, syndromes, zero(errors))
    actual_successes = 0

    for i in axes(syndromes, 2)
      args = decoder.sparse_H, decoder.sparse_HT, syndromes[:, i], 100, fill(decoder.per, decoder.n), zeros(decoder.s, decoder.n), zeros(decoder.s, decoder.n), zeros(decoder.n), zeros(decoder.n)
      decoded_error, converged = LDPCDecoders.syndrome_decode(args...)
      actual_successes += decoded_error == errors[:, i]

      # Test whether the new decoder works similar to the old one
      @assert decoded_error == guess_errors[:, i]
    end

    ler = (num_trials - sum(success)) / num_trials
    return ler
  end

  @test test_bp_decoder()

  # There is a low possibility of these tests failing
  @test test_bp_decoder_batch() < 0.005
  @test test_deprecated_syndrome_decoder() < 0.005
end
