using LDPCDecoders
using JET
using Test

import LinearAlgebra, DelimitedFiles

@testset "JET checks" begin
    rep = report_package("LDPCDecoders";
        ignored_modules=(
            AnyFrameModule(LinearAlgebra),
            AnyFrameModule(DelimitedFiles),
            AnyFrameModule(Base.Broadcast),
        )
    )
    @show rep
    #@test_broken length(JET.get_reports(rep)) == 0
    @test length(JET.get_reports(rep)) == 0
end
