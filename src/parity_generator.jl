using Random

function parity_check_matrix(n, rw, cw)
  n_equations = (n  * cw) ÷ rw
  
  block_size = n_equations ÷ rw
  block = zeros(Int64, block_size, n)
  
  
  for i in 1:block_size
    for j in ((i-1)*cw + 1):((i)*cw)
      block[i,j] = 1
    end
  end
  H = block

  for i in 1:rw - 1
    H = [H; block[:, shuffle(1:end)]] 
  end
  return H
end