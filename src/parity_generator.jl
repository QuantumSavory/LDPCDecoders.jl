"""
parity_check_matrix(n, wr, wc)

Computes a regular LDPC matrix using Callager's algorithm

Block is created initially and then the columns are randomly
permuted to create additional parity check n_equations

# Arguments
* `n`: The length of the code
* `wr`: Row weight (Number of bits a parity check equations acts upon). Must divide n
* `wc`: Column weight (Number of parity check equations for a given bit)

# Examples
```julia
julia> using LDPC
julia> H = parity_check_matrix(1000, 10, 9)
```
"""
function parity_check_matrix(n::Int, wr::Int, wc::Int)

  # For a regular LDPC matrix
  ## wr = wc * (n / n-k)
  @assert n % wr == 0

  n_equations = (n  * wc) ÷ wr
  block_size = n_equations ÷ wc

  block = zeros(Bool, block_size, n)

  for i in 1:block_size
    for j in ((i-1)*wr + 1):((i)*wr)
      block[i,j] = 1
    end
  end

  H = block

  for i in 1:wc - 1
    H = [H; block[:, shuffle(1:end)]]
  end

  return H
end

function save_pcm(H, file_path)
  writedlm(file_path, Int.(H))
end

function load_pcm(file_path)
  H = readdlm(file_path)
  return Int.(H)
end
