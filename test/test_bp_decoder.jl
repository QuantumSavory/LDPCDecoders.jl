@testset "test_bp_decoder.jl" begin
  
  """Test for belief propagation decoder""" function test_bp_decoder()
    H::BitArray{2} = parity_check_matrix(1000, 10, 9)
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
    H::BitArray{2} = parity_check_matrix(1000, 10, 9)
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

    @info "Batch decode logical error rate: $(ler)"
    
    return true
  end

  @test test_bp_decoder()
  @test test_bp_decoder_batch()
end