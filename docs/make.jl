using Documenter
using LDPCDecoders

DocMeta.setdocmeta!(LDPCDecoders, :DocTestSetup, :(using LDPCDecoders, StableRNGs); recursive=true)

makedocs(
    sitename = "LDPCDecoders.jl",
    modules = [LDPCDecoders],
    authors = "QuantumSavory contributors",
    linkcheck = true,
    warnonly = [:missing_docs, :linkcheck],
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
