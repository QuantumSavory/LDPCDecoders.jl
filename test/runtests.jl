using LDPCDecoders
using Test

@testset "LDPCDecoders.jl" begin
    
    ## Syndrome decoding tests
    error_rate = 0.1
    
    pcm = [1 0 1 0 1 0 1; 0 1 1 0 0 1 1; 0 0 0 1 1 1 1]
    error = [0 0 0 0 0 1 0]
    syn = (pcm * error') .% 2
    decoded_error, success = syndrome_decode(pcm, syn, error_rate, 100)
    decoded_error_it, success = syndrome_it_decode(pcm, syn, 100)
    @test decoded_error == vec(error)
    @test decoded_error_it == vec(error)

    pcm = [1 0 1 1 1 0 0; 1 1 0 1 0 1 0; 1 1 1 0 0 0 1]
    error = [0 0 0 1 0 0 0]
    syn = (pcm * error') .% 2
    decoded_error, success = syndrome_decode(pcm, syn, error_rate, 100)
    decoded_error_it, success = syndrome_it_decode(pcm, syn, 100)
    @test decoded_error == vec(error)
    @test decoded_error_it == vec(error)
    
    pcm = [1 1 1 1 0 0 0 0 0 0; 1 0 0 0 1 1 1 0 0 0; 0 1 0 0 1 0 0 1 1 0; 0 0 1 0 0 1 0 1 0 1; 0 0 0 1 0 0 1 0 1 1]
    error = [0 0 0 0 0 0 0 0 0 1] # last bit is an error
    syn = (pcm * error') .% 2
    decoded_error, success = syndrome_decode(pcm, syn, error_rate, 100)
    decoded_error_it, success = syndrome_it_decode(pcm, syn, 100)
    @test decoded_error == vec(error)
    @test decoded_error_it == vec(error)
    

    ## Parity check matrix generation tests
    wr = 10
    wc = 9
    n = 1000
    H = parity_check_matrix(n, wr, wc)
    rsums = sum(H, dims=2)
    csums = sum(H, dims=1)

    @test rsums[1] == wr
    @test csums[1] == wc
    @test all(rsums[1] .== rsums)
    @test all(csums[1] .== csums)
        

end
