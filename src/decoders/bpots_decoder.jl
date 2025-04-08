using SparseArrays


# STATE FOR STORING BP-OTS COMPUTATIONS
#=
The mutable struct stores the computational state during decoding.    
=#
mutable struct BPOTSState
    messages_vc::Dict{Tuple{Int,Int}, Float64}
    messages_cv::Dict{Tuple{Int,Int}, Float64}
    oscillations::Vector{Int}
    biased_nodes::Set{Int}
    prior_decisions::Vector{Int}
    prior_llrs::Vector{Float64}
    
    # Pre-allocated buffers
    Π::Vector{Float64}  # Prior beliefs
    Ω::Vector{Float64}  # Biased beliefs
    decisions::Vector{Int}  # Current decisions
    llrs::Vector{Float64}  # Current LLRs
    best_decisions::Vector{Int}  # Best solution so far
    check_result::Vector{Bool}  # For syndrome checking
    int_check_result::Vector{Int}  # Integer buffer for matrix multiplication
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

    var_neighbors::Vector{Vector{Int}}  # Pre-computed variable neighbors
    check_neighbors::Vector{Vector{Int}}  # Pre-computed check neighbors
end



# BP-OTS STATE INITIALIZATION
# This creates the inital computational state for the BP-OTS algorithm
function initialize_bpots_state(H::SparseMatrixCSC, n::Int)
    # Initializing 2 variables to store messages that flow through the Tanner graph
    messages_vc = Dict{Tuple{Int,Int}, Float64}()
    messages_cv = Dict{Tuple{Int,Int}, Float64}()
    
    # Initialize messages with small random values
    rows = rowvals(H)
    vals = nonzeros(H)
    
    for j in 1:size(H,2)
        for idx in nzrange(H, j)
            i = rows[idx]
            if vals[idx]
                messages_vc[(j,i)] = 0.0
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
        zeros(Float64, n), # prior_llrs
        zeros(Float64, n), # Π
        zeros(Float64, n), # Ω
        zeros(Int, n),     # decisions
        zeros(Float64, n), # llrs
        zeros(Int, n),     # best_decisions
        zeros(Bool, size(H, 1)), # check_result
        zeros(Int, size(H, 1))   # int_check_result (new buffer)
    )
end



# CONSTRUCTOR FOR THE DECODER
# In the BPOTSDecoder constructor, add:
function BPOTSDecoder(H::Union{SparseMatrixCSC{Bool,Int}, BitMatrix}, per::Float64, max_iters::Int; T::Int=9, C::Float64=2.0)
    s, n = size(H)
    sparse_H = sparse(H)
    sparse_HT = sparse(H')
    scratch = initialize_bpots_state(sparse_H, n)
    
    # Pre-compute neighbors for each variable and check node
    var_neighbors = [Vector{Int}() for _ in 1:n]
    check_neighbors = [Vector{Int}() for _ in 1:s]
    
    for j in 1:n
        for idx in nzrange(sparse_H, j)
            i = rowvals(sparse_H)[idx]
            if nonzeros(sparse_H)[idx]
                push!(var_neighbors[j], i)
                push!(check_neighbors[i], j)
            end
        end
    end
    
    return BPOTSDecoder(per, max_iters, s, n, T, C, sparse_H, sparse_HT, scratch, var_neighbors, check_neighbors)
end


# Initialize beliefs with weaker bias
function initialize_beliefs(n::Int, per::Float64)
    # Convert error probability to LLR with weak bias toward no errors
    Π = fill(0.5 * log((1-per)/per), n)  # Weak initial bias
    return Π
end



# COMPUTES BELIEFS AND MAKE DECISIONS
# After msg passing iterations, this function takes all the informations and decides which qubits most likely have the errors
function compute_beliefs!(decoder::BPOTSDecoder, state::BPOTSState, Ω::Vector{Float64})
    fill!(state.decisions, 0)
    fill!(state.llrs, 0)
    
    for j in 1:decoder.n
        # Sum all incoming messages for this variable node
        llr = Ω[j]  # Start with prior belief or bias
        for i in decoder.var_neighbors[j]
            if haskey(state.messages_cv, (i,j))
                llr += state.messages_cv[(i,j)]
            end
        end
        state.llrs[j] = llr
        
        # Decision threshold at 0 - flip bit if LLR is negative
        state.decisions[j] = llr < 0.0 ? 1 : 0
    end
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



# # UPDATE VARIABLE TO CHECK (Eq 1)   [ν^(ℓ)j→i = Ω_j + ∑(i'∈M(v_j)\{c_i}) μ^(ℓ)_i'→j]
# # Variable node j collects evidence from all sources except check node i, and sends this aggregate belief to check node i
function update_variable_to_check!(decoder::BPOTSDecoder, state::BPOTSState, j::Int, i::Int, H::SparseMatrixCSC, Ω::Vector{Float64})
    # Sum all incoming messages except from target check
    msg_sum = 0.0
    for check in decoder.var_neighbors[j]
        if check != i  # Skip the target check node i
            msg = get(state.messages_cv, (check,j), 0.0)
            msg_sum += msg
        end
    end
    
    # Add prior and store
    msg = Ω[j] + msg_sum
    state.messages_vc[(j,i)] = msg
end



# # UPDATE CHECK TO VARIABLE MSG (Eq 2)  [μ^(ℓ)i→j = (-1)^s_i · 2 tanh^(-1)(∏(j'∈N(c_i)\{v_j}) tanh(1/2 ν^(ℓ-1)_j'→i))]
# # Check node i evaluate the parity constraint based on all other connected variables and informs the variable j about what this implies for j's error status
function update_check_to_variable!(decoder::BPOTSDecoder, state::BPOTSState, i::Int, j::Int, H::SparseMatrixCSC, syndrome::Vector{Bool})
    # Compute product of tanh values from other variables
    prod_tanh = 1.0
    MAX_TANH = 0.99999  # For numerical stability
    
    for var in decoder.check_neighbors[i]
        if var != j  # Skip the target variable node j
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
    
    # Compute message using atanh function - SAFETY CHECKS
    if abs(prod_tanh) >= MAX_TANH
        prod_tanh = prod_tanh > 0 ? MAX_TANH : -MAX_TANH
    end
    msg = 2.0 * atanh(prod_tanh)

    # Limit extreme values to prevent overflow
    MAX_MSG = 100.0
    msg = min(MAX_MSG, max(-MAX_MSG, msg))

    state.messages_cv[(i,j)] = msg
end



# DECODE (MAIN ALGORITHM)
function decode!(decoder::BPOTSDecoder, syndrome::Vector{Bool})
    state = decoder.scratch
    reset!(decoder)
    
    # Uses pre-allocated arrays
    state.Π .= log.((1 .- (2*decoder.per/3)) ./ (2*decoder.per/3))
    state.Ω .= state.Π
    
    # Initialize tracking variables
    fill!(state.best_decisions, 0)
    best_mismatch = length(syndrome)
    best_weight = decoder.n
    
    for iter in 1:decoder.max_iters
        # Message passing
        for j in 1:decoder.n
            for i in decoder.var_neighbors[j]
                update_variable_to_check!(decoder, state, j, i, decoder.sparse_H, state.Ω)
            end
        end
        
        for i in 1:decoder.s
            for j in decoder.check_neighbors[i]
                update_check_to_variable!(decoder, state, i, j, decoder.sparse_H, syndrome)
            end
        end
        
        # Compute beliefs in-place
        compute_beliefs!(decoder, state, state.Ω)
        
        # Update oscillations using pre-allocated arrays
        if iter > 1
            for j in 1:decoder.n
                state.oscillations[j] += (state.decisions[j] ⊻ state.prior_decisions[j])
            end
        end
        state.prior_decisions .= state.decisions
        state.prior_llrs .= state.llrs
        
        # Check solution in-place - FIXED to avoid Bool(2) error
        fill!(state.int_check_result, 0)  # Clear the integer buffer
        mul!(state.int_check_result, decoder.sparse_H, state.decisions)  # Matrix multiply into integer buffer
        
        # Safely convert to Bool after applying modulo 2
        mismatch = 0
        for i in 1:length(state.int_check_result)
            # Calculate mod 2 and compare with syndrome
            current_mismatch = (mod(state.int_check_result[i], 2) != syndrome[i])
            state.check_result[i] = current_mismatch
            if current_mismatch
                mismatch += 1
            end
        end
        
        weight = sum(state.decisions)
        
        # Update best solution
        if mismatch < best_mismatch || (mismatch == best_mismatch && weight < best_weight)
            best_mismatch = mismatch
            best_weight = weight
            state.best_decisions .= state.decisions
            
            if mismatch == 0
                return state.best_decisions, true
            end
        end
        
        # Apply biasing if needed - following Algorithm 1 in the paper
        if mismatch > 0 && iter % decoder.T == 0
            # Reset priors to original values
            state.Ω .= state.Π
            
            # Check if there are any oscillating nodes
            if maximum(state.oscillations) > 0
                # Find node with maximum oscillations and minimum LLR
                max_osc = 0
                j1 = 0
                min_llr = Inf
                for j in 1:decoder.n
                    if state.oscillations[j] > max_osc
                        max_osc = state.oscillations[j]
                        j1 = j
                        min_llr = abs(state.llrs[j])
                    elseif state.oscillations[j] == max_osc && abs(state.llrs[j]) < min_llr
                        j1 = j
                        min_llr = abs(state.llrs[j])
                    end
                end
                
                # Reset oscillation counter for j1
                if j1 > 0
                    state.oscillations[j1] = 0
                    
                    # Bias j1
                    state.Ω[j1] = -decoder.C
                end
                
                # Find j2 with minimum |LLR| overall
                j2 = 1
                min_llr = abs(state.llrs[1])
                for j in 2:decoder.n
                    if abs(state.llrs[j]) < min_llr
                        j2 = j
                        min_llr = abs(state.llrs[j])
                    end
                end
                
                # Bias j2 (which may coincide with j1)
                state.Ω[j2] = -decoder.C
            end
        end
    end
    
    return state.best_decisions, false
end


#=

Added pre-allocated buffers to the BPOTSState struct instead of creating new arrays in each iteration
    Added Π and Ω vectors for beliefs (previously created in each decode! call)
    Added decisions and llrs vectors (previously created in each iteration)
    Added best_decisions vector (previously created and copied)
    Added check_result and int_check_result for syndrome checking

Added Pre-computed Neighbor Lists
    var_neighbors: Pre-computed lists of check nodes connected to each variable
    check_neighbors: Pre-computed lists of variable nodes connected to each check
This eliminated repeated calls to get_variable_neighbors() and get_check_neighbors(), which were creating new arrays on each call.

Used .= for assignments instead of = or copy()
Used fill!() to reset arrays rather than creating new ones
Used mul!() for matrix multiplication instead of the * operator  
    (Added an integer buffer for intermediate results, Applied modulo 2 operation explicitly to avoid Boolean conversion errors)
    Since when you multiply a Boolean sparse matrix (decoder.sparse_H) with an integer vector (state.decisions), the result can be greater than 1, and trying bool(2) gives an error.

Replaced findall() and argmin() with direct iteration ()

Modified compute_beliefs! to update pre-allocated buffers

=#