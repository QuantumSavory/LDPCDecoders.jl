using Documenter
using LDPCDecoders

DocMeta.setdocmeta!(LDPCDecoders, :DocTestSetup, :(using LDPCDecoders); recursive=true)
doctest(LDPCDecoders)