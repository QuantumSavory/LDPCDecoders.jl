using LinearAlgebra
using RowEchelon
include("util.jl")

# function parity_to_generator(H::Matrix{Int})
#     rank, n = size(H)
#     n = size(H, 2)
#     At = H[:, 1:n-rank]
#     A = transpose(At)
#     newI  = Matrix{Int}(I, n - rank, n - rank)
#     G = [newI A]
#     return G
# end


function gaussjordan(X)
    """
    Compute binary row reduced echelon
    """

    m, n = size(X)
    P = Matrix(I, m, m)
    
    pivot_old = 0
    for j in 1:n
        filter_down = X[pivot_old+1:m, j]
        pivot = argmax(filter_down) + pivot_old
        println("Argmax ", argmax(filter_down))
        

        if pivot <= m && X[pivot, j] == 1
            pivot_old += 1
            if pivot_old != pivot
                aux = X[pivot, :]
                X[pivot, :] = X[pivot_old, :]
                X[pivot_old, :] = aux 
                
                temp = P[pivot_old, :]
                P[pivot, :] = P[pivot_old, :]
                P[pivot_old, :] = temp
            end

            for i in 1:m
                if i != pivot_old && X[i, j] == 1
                    P[i, :] = abs.(P[i, :] - P[pivot_old, :])
                    X[i, :] = abs.(X[i, :] - X[pivot_old, :])
                end
            end
        end
        
        if pivot_old == m
            break
        end
    end

    return X, P
end

function parity_to_generator(H)
    n_equations, n_code = size(H)


    Href, tQ = gaussjordan(H')
    display(Href)
    Href_diag, aux = gaussjordan(Href')
    tQ = tQ'
    display(tQ)
    display(Href_diag)
    n_bits = n_code - sum(Href_diag)
    println(n_code, n_bits, n_equations)
    Y = zeros(Int8, n_code, n_bits)
    Y[n_code - n_bits + 1:end, :] = Matrix(I, n_bits, n_bits)
    display(Y)

    G = tQ * Y
    return G

end