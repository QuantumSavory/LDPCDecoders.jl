using LinearAlgebra

function generator(H::Matrix{int})
    Ht = transpose(H)
    rref_Ht = rref(Ht)
    G = transpose(rref_Ht[:, end - (size(Ht, 1) - rank(Ht)) + 1:end])
    return G
end
