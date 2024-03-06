using LDPCDecoders
using Test
using SparseArrays

using LDPCDecoders: syndrome_decode

@testset "old LDPCDecoders tests caried over" begin
    ## Syndrome decoding tests
    error_rate = 0.1

    pcm = [1 0 1 0 1 0 1; 0 1 1 0 0 1 1; 0 0 0 1 1 1 1]
    num_bits = 7
    num_checks = 3
    error = [0 0 0 0 0 1 0]
    syn = (pcm * error') .% 2

    args = (sparse(pcm), sparse(pcm'), syn, 10, fill(error_rate, num_bits), zeros(num_checks, num_bits), zeros(num_checks, num_bits), zeros(num_bits), zeros(num_bits))
    decoded_error, success = syndrome_decode(args...)

    @test decoded_error == vec(error)


    pcm = [1 0 1 1 1 0 0; 1 1 0 1 0 1 0; 1 1 1 0 0 0 1]
    num_bits = 7
    num_checks = 3
    error = [0 0 0 1 0 0 0]
    syn = (pcm * error') .% 2

    args = (sparse(pcm), sparse(pcm'), syn, 10, fill(error_rate, num_bits), zeros(num_checks, num_bits), zeros(num_checks, num_bits), zeros(num_bits), zeros(num_bits))
    decoded_error, success = syndrome_decode(args...)
    @test decoded_error == vec(error)

    pcm = [1 1 1 1 0 0 0 0 0 0; 1 0 0 0 1 1 1 0 0 0; 0 1 0 0 1 0 0 1 1 0; 0 0 1 0 0 1 0 1 0 1; 0 0 0 1 0 0 1 0 1 1]
    num_bits = 10
    num_checks = 5
    error = [0 0 0 0 0 0 0 0 0 1] # last bit is an error
    syn = (pcm * error') .% 2

    args = (sparse(pcm), sparse(pcm'), syn, 10, fill(error_rate, num_bits), zeros(num_checks, num_bits), zeros(num_checks, num_bits), zeros(num_bits), zeros(num_bits))
    decoded_error, success = syndrome_decode(args...)

    @test decoded_error == vec(error)

    ## Parity check matrix generation tests
    wr = 10
    wc = 9
    n = 1000
    H = LDPCDecoders.parity_check_matrix(n, wr, wc)
    rsums = sum(H, dims=2)
    csums = sum(H, dims=1)

    @test rsums[1] == wr
    @test csums[1] == wc
    @test all(rsums[1] .== rsums)
    @test all(csums[1] .== csums)

    #TODO:syndrome_it_decode tests
end
