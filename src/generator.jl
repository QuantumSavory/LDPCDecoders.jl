using LinearAlgebra
using RowEchelon
include("util.jl")

function parity_to_generator(H::Matrix{Int})
    rank, n = size(H)
    n = size(H, 2)
    At = H[:, 1:n-rank]
    A = transpose(At)
    newI  = Matrix{Int}(I, n - rank, n - rank)
    G = [newI A]
    return G
end