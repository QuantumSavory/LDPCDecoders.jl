@testitem "JET checks" tags=[:jet] begin
using LDPCDecoders
using JET
using Test

import LinearAlgebra, DelimitedFiles

    rep = report_package("LDPCDecoders";
        ignored_modules=(
            AnyFrameModule(LinearAlgebra),
            AnyFrameModule(DelimitedFiles),
            AnyFrameModule(Base.Broadcast),
        )
    )
    @show rep
    @test length(JET.get_reports(rep)) <= 5
    @test_broken length(JET.get_reports(rep)) == 0
end
