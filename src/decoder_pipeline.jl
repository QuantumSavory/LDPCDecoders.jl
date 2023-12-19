abstract type AbstractSyndromeDecoder end

"""Calls evaluate_classical_decoder for the X and Z checks independently"""
function evaluate_decoder(d::AbstractSyndromeDecoder, nsamples, init_error, gate_error, syndrome_circuit_func, encoding_circuit_func)
    pre_X = [sHadamard(i) for i in n-k+1:n]
    X_error = evaluate_classical_decoder(d, nsamples, init_error, gate_error, syndrome_circuit_func, encoding_circuit_func, logicalxview, 1, d.k, pre_X)
    Z_error = evaluate_classical_decoder(d, nsamples, init_error, gate_error, syndrome_circuit_func, encoding_circuit_func, logicalzview, d.k + 1, 2 * d.k)
    return (X_error, Z_error)
end

"""Uses a Pauli frame simulation to evaluate each part of the CSS code (X or Z tableau)"""
function evaluate_classical_decoder(d::AbstractSyndromeDecoder, nsamples, init_error, gate_error, syndrome_circuit_func, encoding_circuit_func, logical_view_function, guess_start, guess_stop, pre_circuit = nothing)
    H = d.H
    O = d.faults_matrix
    syndrome_circuit = syndrome_circuit_func(H)

    n = d.n
    s = d.s
    k = d.k

    errors = [PauliError(i, init_error) for i in 1:n];

    md = MixedDestabilizer(H)
    
    full_circuit = []

    logview = logical_view_function(md)
    logcirc, _ = syndrome_circuit_func(logview)

    noisy_syndrome_circuit = add_two_qubit_gate_noise(syndrome_circuit, gate_error);
    
    for gate in logcirc
        type = typeof(gate)
        if type == sMRZ
            push!(syndrome_circuit, sMRZ(gate.qubit+s, gate.bit+s))
        else
            push!(syndrome_circuit, type(gate.q1, gate.q2+s))
        end
    end

    ecirc = encoding_circuit_func(syndrome_circuit)
    if isnothing(pre_circuit)
        full_circuit = vcat(pre_circuits, ecirc, errors, noisy_syndrome_circuit)
    else
        full_circuit = vcat(ecirc, errors, noisy_syndrome_circuit)
    end

    frames = PauliFrame(nframes, n+s+k, s+k)
    pftrajectories(frames, full_circuit)
    syndromes = pfmeasurements(frames)[:, 1:s]
    logical_syndromes = pfmeasurements(frames)[:, s+1: s+k]

    for i in 1:nsamples
        guess = decode(d, syndromes[i])

        result = (O * (guess))[guess_start:guess_stop]
        
        if result == logical_syndromes[i]
            decoded += 1
        end
    end

    return (nsamples - decoded) / nsamples
end


"""Type for evaluating the table decoder"""
struct TableDecoder <: AbstractSyndromeDecoder
    H
    faults_matrix
    n
    s
    k
    lookup_table
    time 
end

function TableDecoder(Hx, Hz)
    c = CSS(Hx, Hz)
    H = parity_checks(c)
    s, n = size(H)
    _, _, r = canonicalize!(Base.copy(H), ranks=true)
    k = n - r
    lookup_table, time, _ = @timed create_lookup_table(H)
    faults_matrix = faults_matrix(H)
    return TableDecoder(H, n, s, k, faults_matrix, lookup_table, time)
end


"""Type for evaluating the belief propagation decoder"""
struct BeliefPropDecoder <: AbstractSyndromeDecoder
    H
    faults_matrix
    n
    s
    k
    log_probabs
    channel_probs
    numchecks_X
    b2c_X
    c2b_X
    numchecks_Z
    b2c_Z
    c2b_Z
    err
end

function BeliefPropDecoder(Hx, Hz)
    c = CSS(Hx, Hz)
    H = parity_checks(c)
    s, n = size(H)
    _, _, r = canonicalize!(Base.copy(H), ranks=true)
    k = n - r
    println(typeof(H))
    faults_matrix = faults_matrix(H)
    log_probabs = zeros(n)
    channel_probs = fill(p_init, n)

    numchecks_X = size(Cx)[1]
    b2c_X = zeros(numchecks_X, n)
    c2b_X = zeros(numchecks_X, n)

    numchecks_Z = size(Cz)[1]
    b2c_Z = zeros(numchecks_Z, n)
    c2b_Z = zeros(numchecks_Z, n)
    err = zeros(n)
    return BeliefPropDecoder(H, faults_matrix, n, s, k, log_probabs, channel_probs, numchecks_X, b2c_X, c2b_X, numchecks_Z, b2c_Z, c2b_Z, err)
end


"""Decoding function for the table decoder"""
function decode(d::TableDecoder, syndrome_sample)
    return get(d.lookup_table, syndrome_sample, nothing)
end

"""Decoding function for the belief propagation decoder"""
function decode(d::BeliefPropDecoder, syndrome_sample)
    row_x = syndrome_sample[1:d.numchecks_X]
    row_z = syndrome_sample[d.numchecks_X+1:d.numchecks_X+d.numchecks_Z]

    KguessX, success = syndrome_decode(sparse(d.Cx), sparse(d.Cx'), d.row_x, d.max_iters, d.channel_probs, d.b2c_X, d.c2b_X, d.log_probabs, Base.copy(d.err))
    KguessZ, success = syndrome_decode(sparse(d.Cz), sparse(d.Cz'), d.row_z, d.max_iters, d.channel_probs, d.b2c_Z, d.c2b_Z, d.log_probabs, Base.copy(d.err))
    guess = vcat(KguessZ, KguessX)
end


# Added for running above methods, should not be included in final code(?), implemented inside of QuantumClifford
"""An arbitrary CSS error correcting code defined by its X and Z checks."""
struct CSS <: AbstractECC
    Hx
    Hz
    """Creates a CSS code using the two provided matrices where Hx contains the X checks and Hz contains the Z checks."""
    function CSS(Hx, Hz)
        n = size(Hx, 2)
        if n != size(Hz, 2) error("When constructing a CSS quantum code, the two classical codes are required to have the same block size") end
        if size(Hx,1)+size(Hz,1) >= n error("When constructing a CSS quantum code, the total number of checks (rows) in the parity checks of the two classical codes have to be lower than the block size (the number of columns).") end
        return new(Hx, Hz)
    end
end

function boolean_tableau(c::CSS)
    Hx_height, Hx_width = size(c.Hx)
    Hz_height, Hz_width = size(x.Hz)
    checks_matrix = falses(Hx_height + Hz_height, Hx_width + Hz_width)
    checks_matrix[1:Hx_height, 1:Hx_width] = c.Hx
    checks_matrix[Hx_height+1:end, Hx_width+1:end] = c.Hz
    return CSS(checks_matrix)
end

"""Returns the stabilizer making up the parity check tableau."""
function parity_checks(c::CSS)
    extended_Hx = Matrix{Bool}(vcat(c.Hx, zeros(size(c.Hz))))
    extended_Hz = Matrix{Bool}(vcat(zeros(size(c.Hx)), c.Hz))
    Stabilizer(fill(0x0, size(c.Hx, 2) + size(c.Hz, 2)), extended_Hx, extended_Hz)
end

"""Returns the block length of the code."""
code_n(c::CSS) = size(c.Hx,2)

"""Returns the depth of the parity check matrix"""
code_m(c::CSS) = size(c.Hx, 1) + size(c.Hz, 1)

"""Returns the number of encoded qubits"""
code_k(c::CSS) = (2 * size(c.Hx,2)) - code_m(c)





# NOT COMPLETE, STILL IN PROGRESS
# This will idealy be a simplified implementation based upon evaluating each of the 2 CSS sub-matrices at a time
function evaluate_classical_decoder(H, nsamples, init_error, gate_error, syndrome_circuit_func, encoding_circuit_func, logical_view_func, decoder_func, pre_circuit = nothing)
    decoded = 0

    H_stab = Stabilizer(fill(0x0, size(Hx, 2)), H, zeros(Bool, size(H)))

    O = faults_matrix(H_stab)
    syndrome_circuit = syndrome_circuit_func(H_stab)
    
    s, n = size(H)
    k = n - s

    errors = [PauliError(i, init_error) for i in 1:n];

    md = MixedDestabilizer(H_stab)
    
    full_circuit = []

    logview = logical_view_func(md)
    logcirc, _ = syndrome_circuit_func(logviews)

    noisy_syndrome_circuit = add_two_qubit_gate_noise(syndrome_circuit, gate_error);
    
    for gate in logcirc
        type = typeof(gate)
        if type == sMRZ
            push!(circuit, sMRZ(gate.qubit+s, gate.bit+s))
        else
            push!(circuit, type(gate.q1, gate.q2+s))
        end
    end

    ecirc = encoding_circuit_func(syndrome_circuit)
    if isnothing(pre_circuit)
        full_circuit = vcat(pre_circuits, ecirc, errors, noisy_syndrome_circuit)
    else
        full_circuit = vcat(ecirc, errors, noisy_syndrome_circuit)
    end

    frames = PauliFrame(nframes, n+s+k, s+k)
    pftrajectories(frames, full_circuit)
    syndromes = pfmeasurements(frames)[:, 1:s]
    logical_syndromes = pfmeasurements(frames)[:, s+1: s+k]

    for i in 1:nsamples
        guess = decode(decoder_obj, syndromes[i]) # TODO: replace 'decoder_obj' with proper object

        result = (O * (guess))
        
        if result == logical_syndromes[i]
            decoded += 1
        end
    end

    return (nsamples - decoded) / nsamples
end