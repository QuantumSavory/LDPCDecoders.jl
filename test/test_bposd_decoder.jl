using Test
using LDPCDecoders

@testset "test_bposd_decoder.jl" begin

  """Test for BP-OSD decoder"""
  function test_bposd_decoder()
    H = LDPCDecoders.parity_check_matrix(1000, 10, 9)
    per = 0.01
    err = rand(1000) .< per
    syn = (H * err) .% 2

    bposd = BeliefPropagationOSDDecoder(H, per, 100)
    guess, success = decode!(bposd, syn)

    return guess == err
  end

  """Test high order OSD"""
  function test_bposd_decoder_high_order()
    H = LDPCDecoders.parity_check_matrix(1000, 10, 9)
    per = 0.01
    err = rand(1000) .< per
    syn = (H * err) .% 2

    orders = 2:5
    succ = true
    for osd_order in orders
      bposd = BeliefPropagationOSDDecoder(H, per, 100; osd_order=osd_order)
      guess, success = decode!(bposd, syn)
      succ = succ & (guess == err)
    end

    return succ
  end

  """Test for BP-OSD decoder with large error rate. Even if the decoding is not accurate, OSD will still ensure consistency between guess and syndromes."""
  function test_bposd_decoder_large_error_rate()
    H = LDPCDecoders.parity_check_matrix(1000, 10, 9)
    per = 0.2
    err = rand(1000) .< per
    syn = (H * err) .% 2

    bposd = BeliefPropagationOSDDecoder(H, per, 100)
    guess, success = decode!(bposd, syn)

    return syn == (H * guess) .% 2
  end

  @test test_bposd_decoder()
  @test test_bposd_decoder_high_order()
  @test test_bposd_decoder_large_error_rate()
end
