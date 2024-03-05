@testset "test_bp_decoder.jl" begin
  
  """Test for belief propagation decoder""" 
  function test_bp_decoder()
    H::BitArray{2} = LDPCDecoders.parity_check_matrix(1000, 10, 9)
    per = 0.01
    err::BitArray{1} = rand(1000) .< per 
    syn::BitArray{1} = (H * err) .% 2
    
    bpd = BeliefPropagationDecoder(H, per, 100)
    decoded_error::BitArray{1} = zeros(1000)
  
    success = decode!(bpd, syn, decoded_error)

    return success
  end

  """Test for batch belief propagation decoder"""
  function test_bp_decoder_batch()
    H::BitArray{2} = LDPCDecoders.parity_check_matrix(1000, 10, 9)
    per = 0.1
    num_trials = 100
    errors::BitArray{2} = rand(Base.Float64, (num_trials,1000)) .< per
    syndromes::BitArray{2} = zeros(num_trials, 900)
    
    for i in axes(errors, 1)
      syn = (H * errors[i, :]) .% 2
      syndromes[i, :] = syn
    end    

    bpd = BeliefPropagationDecoder(H, per, 100)
    success::BitArray{1} = batchdecode!(bpd, syndromes, errors)
    ler = (num_trials - sum(success)) / num_trials
    return ler
  end

  """Test for batch belief propagation decoder"""
  function test_depricated_syndrome_decoder()
    H::BitArray{2} = LDPCDecoders.parity_check_matrix(1000, 10, 9)
    per = 0.1
    num_trials = 100
    errors::BitArray{2} = rand(Base.Float64, (num_trials,1000)) .< per
    syndromes::BitArray{2} = zeros(num_trials, 900)
    
    for i in axes(errors, 1)
      syn = (H * errors[i, :]) .% 2
      syndromes[i, :] = syn
    end    
    
    decoder = BeliefPropagationDecoder(H, per, 100)
    success::BitArray{1} = batchdecode!(decoder, syndromes, errors)
    
    for i in axes(syndromes, 1)
      args = decoder.sparse_H, decoder.sparse_HT, syndromes[i, :], 100, fill(decoder.per, decoder.n), zeros(decoder.s, decoder.n), zeros(decoder.s, decoder.n), zeros(decoder.n), zeros(decoder.n)
      decoded_error, converged = syndrome_decode(args...)
      
      # Test whether the new decoder works similar to the old one
      @assert decoded_error == errors[i, :]
    end

    ler = (num_trials - sum(success)) / num_trials
    return ler
  end

  @test test_bp_decoder()

  # There is a low possibility of these tests failing
  @test test_bp_decoder_batch() < 0.9
  @test test_depricated_syndrome_decoder() < 0.9
end