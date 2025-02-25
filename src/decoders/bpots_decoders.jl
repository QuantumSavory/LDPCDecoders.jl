using SparseArrays

# State for storing BP-OTS computations
mutable struct BPOTSState
    messages_vc::Dict{Tuple{Int,Int}, Float64}
    messages_cv::Dict{Tuple{Int,Int}, Float64}
    oscillations::Vector{Int}
    biased_nodes::Set{Int}
    prior_decisions::Vector{Int}
    prior_llrs::Vector{Float64}
end

# Main decoder struct
struct BPOTSDecoder <: AbstractDecoder
    per::Float64        # Physical error rate
    max_iters::Int      # Maximum iterations
    s::Int             # Number of stabilizers
    n::Int             # Number of qubits
    T::Int             # Biasing period
    C::Float64         # Bias constant
    sparse_H::SparseMatrixCSC{Bool,Int}  # Sparse parity check matrix
    sparse_HT::SparseMatrixCSC{Bool,Int} # Transposed matrix
    scratch::BPOTSState # Working space for computations
end

# BP-OTS State initialization
# Initialize state
function initialize_bpots_state(H::SparseMatrixCSC, n::Int)
    messages_vc = Dict{Tuple{Int,Int}, Float64}()
    messages_cv = Dict{Tuple{Int,Int}, Float64}()
    
    # Initialize messages with small random values
    rows = rowvals(H)
    vals = nonzeros(H)
    
    for j in 1:size(H,2)
        for idx in nzrange(H, j)
            i = rows[idx]
            if vals[idx]
                messages_vc[(j,i)] = 0.0  # Start neutral
                messages_cv[(i,j)] = 0.0
            end
        end
    end
    
    return BPOTSState(
        messages_vc,
        messages_cv,
        zeros(Int, n),     # oscillations
        Set{Int}(),        # biased_nodes
        zeros(Int, n),     # prior_decisions
        zeros(Float64, n)  # prior_llrs
    )
end

# Helper to get check node neighbors
function get_check_neighbors(H::SparseMatrixCSC, check_idx::Int)
    neighbors = Int[]
    # Look through the check node's row
    for col in 1:size(H,2)
        if H[check_idx, col]
            push!(neighbors, col)
        end
    end
    return neighbors
end

# Helper to get variable node neighbors 
function get_variable_neighbors(H::SparseMatrixCSC, var_idx::Int)
    neighbors = Int[]
    # Look through column's non-zero entries
    for idx in nzrange(H, var_idx)
        row = rowvals(H)[idx]
        if nonzeros(H)[idx]
            push!(neighbors, row)
        end
    end
    return neighbors
end

# Constructor for the decoder
function BPOTSDecoder(H::Union{SparseMatrixCSC{Bool,Int}, BitMatrix}, per::Float64, max_iters::Int; T::Int=9, C::Float64=2.0)
    s, n = size(H)
    sparse_H = sparse(H)
    sparse_HT = sparse(H')
    scratch = initialize_bpots_state(sparse_H, n)
    
    return BPOTSDecoder(per, max_iters, s, n, T, C, sparse_H, sparse_HT, scratch)
end

# Initialize beliefs with weaker bias
function initialize_beliefs(n::Int, per::Float64)
    # Convert error probability to LLR with weak bias toward no errors
    Π = fill(0.5 * log((1-per)/per), n)  # Weak initial bias
    return Π
end

# Compute beliefs and make decisions
function compute_beliefs!(decoder::BPOTSDecoder, state::BPOTSState, Ω::Vector{Float64})
    n = decoder.n
    decisions = zeros(Int, n)
    llrs = zeros(Float64, n)
    
    for j in 1:n
        # Sum all incoming messages for this variable node
        llr = Ω[j]  # Start with prior
        for i in get_variable_neighbors(decoder.sparse_H, j)
            if haskey(state.messages_cv, (i,j))
                llr += state.messages_cv[(i,j)]
            end
        end
        llrs[j] = llr
        
        # Decision threshold at 0 - flip bit if LLR is negative
        decisions[j] = llr < 0.0 ? 1 : 0
    end
    
    return decisions, llrs
end

# Reset state between decodings with proper initialization
function reset!(decoder::BPOTSDecoder)
    state = decoder.scratch
    empty!(state.biased_nodes)
    fill!(state.oscillations, 0)
    fill!(state.prior_decisions, 0)
    fill!(state.prior_llrs, 0)
    
    # Reset messages to neutral values (paper starts with 0)
    for k in keys(state.messages_vc)
        state.messages_vc[k] = 0.0
    end
    for k in keys(state.messages_cv)
        state.messages_cv[k] = 0.0
    end
    
    return decoder
end

# Update variable-to-check message according to Equation (1) in the paper
function update_variable_to_check!(state::BPOTSState, j::Int, i::Int, H::SparseMatrixCSC, Ω::Vector{Float64})
    connected_checks = get_variable_neighbors(H, j)
    
    # Sum all incoming messages except from target check
    msg_sum = 0.0
    for check in connected_checks
        if check != i
            msg = get(state.messages_cv, (check,j), 0.0)
            msg_sum += msg
        end
    end
    
    # Add prior and store - Equation (1) from the paper
    msg = Ω[j] + msg_sum
    state.messages_vc[(j,i)] = msg
end

# Update check-to-variable message according to Equation (2) in the paper
function update_check_to_variable!(state::BPOTSState, i::Int, j::Int, H::SparseMatrixCSC, syndrome::Vector{Bool})
    connected_vars = get_check_neighbors(H, i)
    
    # Compute product of tanh values from other variables
    prod_tanh = 1.0
    MAX_TANH = 0.99999  # For numerical stability
    
    for var in connected_vars
        if var != j
            msg = get(state.messages_vc, (var,i), 0.0)
            t = tanh(0.5 * msg)
            t = min(MAX_TANH, max(-MAX_TANH, t))
            prod_tanh *= t
        end
    end
    
    # Apply syndrome as in Equation (2)
    if syndrome[i]
        prod_tanh = -prod_tanh
    end
    
    # Compute message using atanh function
    msg = 2.0 * atanh(prod_tanh)
    state.messages_cv[(i,j)] = msg
end

# decode! function matching the paper's algorithm
function decode!(decoder::BPOTSDecoder, syndrome::Vector{Bool})
    state = decoder.scratch
    reset!(decoder)
    
    println("\nDEBUG: Starting decode")
    println("DEBUG: Syndrome: ", syndrome)
    
    # Initialize priors
    # Πj = log((1-(2ϵ/3))/(2ϵ/3))- from the paper
    Π = fill(log((1-(2*decoder.per/3))/(2*decoder.per/3)), decoder.n)
    Ω = copy(Π)
    
    # Track best solution
    best_decisions = zeros(Int, decoder.n)
    best_mismatch = length(syndrome)
    best_weight = decoder.n
    
    for iter in 1:decoder.max_iters
        println("\nDEBUG: === Iteration $iter ===")
        
        # Message passing
        println("\nDEBUG: Variable updates")
        for j in 1:decoder.n
            for i in get_variable_neighbors(decoder.sparse_H, j)
                update_variable_to_check!(state, j, i, decoder.sparse_H, Ω)
            end
        end
        
        println("\nDEBUG: Check updates")
        for i in 1:decoder.s
            for j in 1:decoder.n
                if decoder.sparse_H[i,j]
                    update_check_to_variable!(state, i, j, decoder.sparse_H, syndrome)
                end
            end
        end
        
        # Compute marginals (beliefs)
        decisions = zeros(Int, decoder.n)
        llrs = zeros(Float64, decoder.n)
        
        println("\nDEBUG: Computing beliefs")
        for j in 1:decoder.n
            llr = Ω[j]
            println("DEBUG: Node $j")
            println("DEBUG:   Prior: $llr")
            
            for i in get_variable_neighbors(decoder.sparse_H, j)
                msg = state.messages_cv[(i,j)]
                llr += msg
                println("DEBUG:   Added check $i msg: $msg")
            end
            
            llrs[j] = llr
            decisions[j] = llr < 0.0 ? 1 : 0
            println("DEBUG:   Final LLR: $llr -> Decision: $(decisions[j])")
        end
        
        # Update oscillation vector using XOR between consecutive decisions
        # This tracks how often bits flip
        if iter > 1
            for j in 1:decoder.n
                state.oscillations[j] += (decisions[j] ⊻ state.prior_decisions[j])
            end
        end
        state.prior_decisions = copy(decisions)
        state.prior_llrs = copy(llrs)
        
        # Check solution
        check_result = Bool.(mod.(decoder.sparse_H * decisions, 2))
        mismatch = count(check_result .!= syndrome)
        weight = sum(decisions)
        println("\nDEBUG: Mismatch: $mismatch, Weight: $weight")
        
        if mismatch < best_mismatch || (mismatch == best_mismatch && weight < best_weight)
            println("DEBUG: New best solution!")
            best_mismatch = mismatch
            best_weight = weight
            best_decisions = copy(decisions)
            
            if mismatch == 0
                println("DEBUG: Found valid solution!")
                return best_decisions, true
            end
        end
        
        # Apply biasing if needed - following Algorithm 1 in the paper
        if mismatch > 0 && iter % decoder.T == 0
            # Reset priors to original values
            Ω = copy(Π)
            
            # Check if there are any oscillating nodes
            if maximum(state.oscillations) > 0
                # Step 1: Find nodes with maximum oscillations (F set in the paper)
                max_osc = maximum(state.oscillations)
                max_osc_indices = findall(o -> o == max_osc, state.oscillations)
                
                # Step 2: From F, select j1 with minimum |LLR|
                j1 = max_osc_indices[argmin([abs(llrs[j]) for j in max_osc_indices])]
                
                # Reset oscillation counter for j1
                state.oscillations[j1] = 0
                
                # Bias j1 
                Ω[j1] = -decoder.C
                println("DEBUG: Biasing j1 (node $j1): $(Ω[j1])")
                
                # Step 3: Select j2 with minimum |LLR| overall
                j2 = argmin(abs.(llrs))
                
                # Bias j2 (which may coincide with j1)
                Ω[j2] = -decoder.C
                println("DEBUG: Biasing j2 (node $j2): $(Ω[j2])")
            end
        end
    end
    
    return best_decisions, false
end
