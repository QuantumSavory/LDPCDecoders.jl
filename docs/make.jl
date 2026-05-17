using Documenter
using LDPCDecoders

DocMeta.setdocmeta!(LDPCDecoders, :DocTestSetup, :(using LDPCDecoders); recursive=true)

makedocs(
    sitename = "LDPCDecoders.jl",
    modules = [LDPCDecoders],
    warnonly = [:missing_docs],
    authors = "QuantumSavory contributors",
    pages = [
        "Home" => "index.md",
        "Decoders" => [
            "Belief Propagation" => "decoders/belief_propagation.md",
            "Belief Propagation + OSD" => "decoders/bp_osd.md",
            "Iterative BitFlip" => "decoders/bitflip.md",
            "BP-OTS" => "decoders/bpots.md",
        ],
        "API Reference" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/QuantumSavory/LDPCDecoders.jl.git",
)
