struct BeliefPropagationOSDDecoder <: AbstractDecoder
    """A belief propagation decoder as a subroutine"""
    bp_decoder::BeliefPropagationDecoder
    """Dense form of the parity check matrix"""
    H::BitMatrix
    """The order of OSD; defaulted to be 0 in the constructor"""
    osd_order::Int
end

function BeliefPropagationOSDDecoder(H, per::Float64, max_iters::Int; osd_order::Int=0)
    bp_decoder = BeliefPropagationDecoder(H, per, max_iters)
    return BeliefPropagationOSDDecoder(bp_decoder, H, osd_order)
end

function rowswap!(H::BitMatrix, i, j)
    @inbounds H[i, :], H[j, :] = H[j, :], H[i, :] # TODO This could be further optimized?
end

function decode!(decoder::BeliefPropagationOSDDecoder, syndrome::AbstractVector)
    # use BP to get hard and soft decisions
    bp_err, converged = decode!(decoder.bp_decoder, syndrome) # hard decisions
    bp_log_probabs = decoder.bp_decoder.scratch.log_probabs # soft decisions
    bp_probabs = exp.(bp_log_probabs)
    # sort columns by reliability, less reliable columns first
    sort_by_reliability = sortperm(max.(bp_probabs, 1 .- bp_probabs), rev=true)
    H_sorted = decoder.H[:, sort_by_reliability]
    bp_err_sorted = bp_err[sort_by_reliability]
    # TODO an optimized version of OSD can be implemented when osd_order = 0, see Algorithm 2 in	https://doi.org/10.22331/q-2021-11-22-585
    err = osd(H_sorted, syndrome, bp_err_sorted, decoder.osd_order)
    return err[invperm(sort_by_reliability)], converged # also return whether BP is converged
end

function osd(H, syndrome, bp_err, osd_order)
    m, n = size(H)
    # diagnolize the submatrix corresponding to independent columns via Gaussian elimination
    # first obtain the row canonical form
    # and find least reliable indices, i.e., the first r pivot columns (assume H is rearranged by reliability)
    least_reliable_rows = [] # row indices of pivot elements
    least_reliable_cols = [] # column indices of pivot elements
    r = 0 # compute rank of H
    i, j = 1, 1
    s = copy(syndrome) # transform syndrome along with H in Gaussian elimination

    while i <= m && j <= n
        k = findfirst(H[i:end, j])
        if isnothing(k) # not an independent column
            j += 1
        else
            if k > 1
                ii = i + k - 1 # the first row after `i` with 1 in column `j`
                rowswap!(H, i, ii) # TODO For optimization: Is this swap necessary? We may just track the row index
                s[i], s[ii] = s[ii], s[i]
            end
            for ii in i+1:m
                if H[ii, j]
                    H[ii, :] .⊻= H[i, :]
                    s[ii] ⊻= s[i]
                end
            end
            push!(least_reliable_rows, i)
            push!(least_reliable_cols, j)
            i += 1
            j += 1
            r += 1
        end
    end

    # then obtain a diagonal submatrix on the least reliable part
    for (i, j) in zip(reverse(least_reliable_rows), reverse(least_reliable_cols))
        for ii in 1:i-1
            if H[ii, j]
                H[ii, :] .⊻= H[i, :]
                s[ii] ⊻= s[i]
            end
        end
    end

    if osd_order > n - r
        @warn "The order of OSD $osd_order is greater than the size of the information set $(n-r). We set osd_order = $(n-r)."
        osd_order = n - r
    end

    best_err = copy(bp_err)
    err = Bool.(copy(bp_err)) # TODO why error is in Float in BP?
    most_reliable_cols = setdiff(1:n, least_reliable_cols)
    min_weight = n + 1

    for x in 0:2^osd_order-1
        # first compute the `most_reliable_cols` part of errors
        # try all possible errors on the first `osd_order` bits within `most_reliable_cols`
        if x != 0
            trial_err = BitArray([x >> i & 1 for i in 0:osd_order-1])
            err[most_reliable_cols[1:osd_order]] = trial_err
        end
        # then based on the `most_reliable_cols` part of errors, compute the `least_reliable_cols` part of errors
        for (i, j) in zip(least_reliable_rows, least_reliable_cols)
            err[j] = s[i]
            for k in most_reliable_cols
                err[j] ⊻= H[i, k] * err[k]
            end
        end
        weight = sum(err) # This weight is set for depolarizing noise
        # TODO More generally, it should be a function depending on the noise model
        if weight < min_weight
            min_weight = weight
            best_err = copy(err)
        end
    end

    return best_err
end
