@testitem "Parity check matrix generation" begin
using LDPCDecoders
using Test

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
end
