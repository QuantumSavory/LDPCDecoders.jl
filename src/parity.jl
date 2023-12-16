function hamming_to_parity(rank)
    num_rows = 2^rank - 1

    parity = zeros(Int, num_rows, rank)

    for i in range(1, num_rows)
        digits = base2(i, rank)

        for j in range(1, rank)
            parity[i, j] = digits[j]
        end
    end

    return transpose(parity)
end

function repetition_to_parity(distance)
    parity = zeros(distance - 1, distance)

    for i in range(1, distance - 1)
        parity[i, i] = 1
        parity[i, i + 1] = 1
    end

    return parity
end
