using SparseArrays


# STATE FOR STORING BP-OTS COMPUTATIONS
#=
The mutable struct stores the computational state during decoding.    
=#
mutable struct BPOTSState
    messages_vc::Dict{Tuple{Int,Int}, Float64}  # dictionary mapping (variable, check) pairs to message values
    messages_cv::Dict{Tuple{Int,Int}, Float64} # mapping (check, variable) pairs to message values
    oscillations::Vector{Int} # this is a counter that tracks how many times each bit's decision changes
    biased_nodes::Set{Int} # variable nodes that have been biased
    prior_decisions::Vector{Int} # bit decisions from the previous Iteration
    prior_llrs::Vector{Float64} # log likelihood ratios from the previous iteration
end



# MAIN DECODER STRUCT
struct BPOTSDecoder <: AbstractDecoder
    per::Float64        # Physical error rate (probability of a qubit error)
    max_iters::Int      # Maximum # of iterations before giving up
    s::Int             # Number of stabilizers  (Dimension of the parity check matrix- s x n)
    n::Int             # Number of qubits
    T::Int             # Biasing period (How often to apply the bias)
    C::Float64         # Bias constant
    sparse_H::SparseMatrixCSC{Bool,Int}  # Sparse parity check matrix
    sparse_HT::SparseMatrixCSC{Bool,Int} # Transposed parity check matrix
    scratch::BPOTSState # Working state for computations
end



# BP-OTS STATE INITIALIZATION
# This creates the inital computational state for the BP-OTS algorithm
function initialize_bpots_state(H::SparseMatrixCSC, n::Int)

    # Initializing 2 variables to sotre messages that flow through the Tanner graph
    messages_vc = Dict{Tuple{Int,Int}, Float64}()
    messages_cv = Dict{Tuple{Int,Int}, Float64}()
    
    # Initialize messages with small random values
    rows = rowvals(H) # Since H is sparse, this gets the row indices of the non-zero elements
    vals = nonzeros(H) # Column indices of the non-zero elements
    
    # Loops throught every column j, and for each connected row i, initalizes the msg between those nodes to 0.0
    for j in 1:size(H,2)
        for idx in nzrange(H, j)
            i = rows[idx]
            if vals[idx]
                messages_vc[(j,i)] = 0.0 # Makes sure there is no bias
                messages_cv[(i,j)] = 0.0
            end
        end
    end

    # Create the state and returns it
    return BPOTSState(
        messages_vc,
        messages_cv,
        zeros(Int, n),     # oscillations
        Set{Int}(),        # biased_nodes
        zeros(Int, n),     # prior_decisions
        zeros(Float64, n)  # prior_llrs
    )
end



# CONSTRUCTOR FOR THE DECODER
function BPOTSDecoder(H::Union{SparseMatrixCSC{Bool,Int}, BitMatrix}, per::Float64, max_iters::Int; T::Int=9, C::Float64=2.0)
    s, n = size(H) # Assigns the dimensions to s and n
    sparse_H = sparse(H) # Makes H sparse, incase it is Dense
    sparse_HT = sparse(H') # Transposes H
    scratch = initialize_bpots_state(sparse_H, n) # Creates a workspace to carry out the decoding
    
    return BPOTSDecoder(per, max_iters, s, n, T, C, sparse_H, sparse_HT, scratch)
end



# Initialize beliefs with weaker bias
function initialize_beliefs(n::Int, per::Float64)
    # Convert error probability to LLR with weak bias toward no errors
    Π = fill(0.5 * log((1-per)/per), n)  # Weak initial bias
    return Π
end



# COMPUTER BELIEFS AND MAKE DECISIONS
# After msg passing iterations, this function takes all the informations and decides which qubits most likely have the errors
function compute_beliefs!(decoder::BPOTSDecoder, state::BPOTSState, Ω::Vector{Float64})
    n = decoder.n
    decisions = zeros(Int, n)
    llrs = zeros(Float64, n)
    
    for j in 1:n
        # Sum all incoming messages for this variable node
        llr = Ω[j]  # Start with prior belief or bias
        for i in get_variable_neighbors(decoder.sparse_H, j) # finds all check nodes connected to this variable
            if haskey(state.messages_cv, (i,j))
                llr += state.messages_cv[(i,j)] # adds all the cv msgs to the prior
            end
        end
        llrs[j] = llr
        
        # Decision threshold at 0 - flip bit if LLR is negative
        decisions[j] = llr < 0.0 ? 1 : 0 # 1 if LLR < 0, 0 otherwise
    end
    
    return decisions, llrs
end



#RESEST STATE BETWEEN DECODINGS
#Prepares the decoder for reuase by clearing its internal state
function reset!(decoder::BPOTSDecoder)
    state = decoder.scratch # reference to scratch state object
    empty!(state.biased_nodes) # clears set of nodes with their biases
    fill!(state.oscillations, 0) # resets bit-flip counter to 0
    fill!(state.prior_decisions, 0) # clears the past decisions
    fill!(state.prior_llrs, 0) # clears previous llrs
    
    # Reset messages to neutral values
    for k in keys(state.messages_vc)
        state.messages_vc[k] = 0.0
    end
    for k in keys(state.messages_cv)
        state.messages_cv[k] = 0.0
    end
    
    return decoder
end



# UPDATE VARIABLE TO CHECK (Eq 1)   [ν^(ℓ)j→i = Ω_j + ∑(i'∈M(v_j)\{c_i}) μ^(ℓ)_i'→j]
# Variable node j collects evidence from all sources except check node i, and sends this aggregate belief to check node i
function update_variable_to_check!(state::BPOTSState, j::Int, i::Int, H::SparseMatrixCSC, Ω::Vector{Float64})
    connected_checks = get_variable_neighbors(H, j) # find all check nodes connected to variable j
    
    # Sum all incoming messages except from target check
    msg_sum = 0.0
    for check in connected_checks
        if check != i # skip the target check node i to prevent echo chammer effects
            msg = get(state.messages_cv, (check,j), 0.0) # get the existing msg from the check node to the variable node
            msg_sum += msg
        end
    end
    
    # Add prior and store - Equation (1) from the paper
    msg = Ω[j] + msg_sum #add prior belief/ bias to the sum
    state.messages_vc[(j,i)] = msg
end



# UPDATE CHECK TO VARIABLE MSG (Eq 2)  [μ^(ℓ)i→j = (-1)^s_i · 2 tanh^(-1)(∏(j'∈N(c_i)\{v_j}) tanh(1/2 ν^(ℓ-1)_j'→i))]
# Check node i evaluate the parity constraint based on all other connected variables and informs the variable j about what this implies for j's error status
function update_check_to_variable!(state::BPOTSState, i::Int, j::Int, H::SparseMatrixCSC, syndrome::Vector{Bool})
    connected_vars = get_check_neighbors(H, i) # finds all variable nodes connected to check i
    
    # Compute product of tanh values from other variables
    prod_tanh = 1.0
    MAX_TANH = 0.99999  # For numerical stability
    
    for var in connected_vars
        if var != j # for each connected variable node, this skips the target variable node j
            msg = get(state.messages_vc, (var,i), 0.0) # gets the msg from that variable to this check
            t = tanh(0.5 * msg)
            t = min(MAX_TANH, max(-MAX_TANH, t))
            prod_tanh *= t
        end
    end
    
    # Apply syndrome as in Equation (2)
    if syndrome[i]
        prod_tanh = -prod_tanh # sign flip if syndrome is true
    end
    
    # Compute message using atanh function - SAFETY CHECKS
    # Avoid atanh on values too close to ±1
    if abs(prod_tanh) >= MAX_TANH
        prod_tanh = prod_tanh > 0 ? MAX_TANH : -MAX_TANH
    end
    msg = 2.0 * atanh(prod_tanh)

    # Limit extreme values to prevent overflow
    MAX_MSG = 100.0
    msg = min(MAX_MSG, max(-MAX_MSG, msg))

    state.messages_cv[(i,j)] = msg
end



#=
HELPER FUNCTIONS FOR MESSAGE PASSING
=#

# HELPER TO GET CHECK NODE NEIGHBORS
# This function identifies all variable nodes (qubits) connected to a specific check node (stabilizer) in the tanner graph
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


# HELPER TO GET VARIABLE NODE NEIGHBORS
# This function finds all check nodes (stabilizers) that involve a specific variable node (qubit)
# Since, for msg passing, we need to to know which stabilizers are affected by each qubit
function get_variable_neighbors(H::SparseMatrixCSC, var_idx::Int)
    neighbors = Int[] # To sdtore connected check nodes
    # Look through column's non-zero entries
    for idx in nzrange(H, var_idx) # where var_idx is the col we are examining
        row = rowvals(H)[idx]
        if nonzeros(H)[idx]
            push!(neighbors, row)
        end
    end
    return neighbors
end

#=
END OF HELPPER FUNCTIONS
=#



# DECODE (MAIN ALGORITHM)
function decode!(decoder::BPOTSDecoder, syndrome::Vector{Bool})

    state = decoder.scratch
    reset!(decoder)
    
    #println("\nDEBUG: Starting decode")
    #println("DEBUG: Syndrome: ", syndrome)
    
    # Initialize priors beliefs based on the physical error rate per
    # Πj = log((1-(2ϵ/3))/(2ϵ/3)) from the paper
    Π = fill(log((1-(2*decoder.per/3))/(2*decoder.per/3)), decoder.n)
    Ω = copy(Π)
    
    # Track best solution
    best_decisions = zeros(Int, decoder.n)
    best_mismatch = length(syndrome)
    best_weight = decoder.n
    
    for iter in 1:decoder.max_iters
        #println("\nDEBUG: === Iteration $iter ===")
        
        # Message passing
        #println("\nDEBUG: Variable updates")
        for j in 1:decoder.n
            for i in get_variable_neighbors(decoder.sparse_H, j)
                update_variable_to_check!(state, j, i, decoder.sparse_H, Ω)
            end
        end
        
        #println("\nDEBUG: Check updates")
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
        
        #println("\nDEBUG: Computing beliefs")
        for j in 1:decoder.n
            llr = Ω[j]
            #println("DEBUG: Node $j")
            #println("DEBUG:   Prior: $llr")
            
            for i in get_variable_neighbors(decoder.sparse_H, j)
                msg = state.messages_cv[(i,j)]
                llr += msg
                #println("DEBUG:   Added check $i msg: $msg")
            end
            
            llrs[j] = llr
            decisions[j] = llr < 0.0 ? 1 : 0
            #println("DEBUG:   Final LLR: $llr -> Decision: $(decisions[j])")
        end
        
        # Update oscillation vector using XOR between consecutive decisions
        # This tracks how often bits flip, as described in the paper
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
        #println("\nDEBUG: Mismatch: $mismatch, Weight: $weight")
        
        if mismatch < best_mismatch || (mismatch == best_mismatch && weight < best_weight)
            #println("DEBUG: New best solution!")
            best_mismatch = mismatch
            best_weight = weight
            best_decisions = copy(decisions)
            
            if mismatch == 0
                #println("DEBUG: Found valid solution!")
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
                #println("DEBUG: Biasing j1 (node $j1): $(Ω[j1])")
                
                # Step 3: Select j2 with minimum |LLR| overall
                j2 = argmin(abs.(llrs))
                
                # Bias j2 (which may coincide with j1)
                Ω[j2] = -decoder.C
                #println("DEBUG: Biasing j2 (node $j2): $(Ω[j2])")
            end
        end
    end
    
    return best_decisions, false
end
