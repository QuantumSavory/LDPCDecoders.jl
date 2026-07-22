@testitem "Doctests" tags=[:doctests] begin
using Documenter
using LDPCDecoders
using StableRNGs

ENV["LINES"] = 80    # for forcing `displaysize(io)` to be big enough
ENV["COLUMNS"] = 80
DocMeta.setdocmeta!(LDPCDecoders, :DocTestSetup, :(using LDPCDecoders;); recursive=true)
doctest(LDPCDecoders)
end
