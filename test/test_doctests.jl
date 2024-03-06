using Documenter
using LDPCDecoders

ENV["LINES"] = 80    # for forcing `displaysize(io)` to be big enough
ENV["COLUMNS"] = 80
@testset "Doctests" begin
    DocMeta.setdocmeta!(LDPCDecoders, :DocTestSetup, :(using LDPCDecoders;); recursive=true)
    doctest(LDPCDecoders)
end
