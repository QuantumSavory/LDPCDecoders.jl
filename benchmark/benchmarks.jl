using BenchmarkTools
using LDPCDecoders

const SUITE = BenchmarkGroup()

# ── Setup ────────────────────────────────────────────────────────────────────

H = LDPCDecoders.parity_check_matrix(1000, 10, 9)
per = 0.01
err = rand(1000) .< per
syn = Bool.((H * err) .% 2)

# ── BP-OSD decode! ───────────────────────────────────────────────────────────

SUITE["bposd"] = BenchmarkGroup(["bposd"])
bposd0 = BeliefPropagationOSDDecoder(BitMatrix(Matrix(H)), per, 100; osd_order=0)
SUITE["bposd"]["decode_osd0"] = @benchmarkable decode!(d, s) setup=(d=deepcopy($bposd0); s=copy($syn)) evals=1
bposd2 = BeliefPropagationOSDDecoder(BitMatrix(Matrix(H)), per, 100; osd_order=2)
SUITE["bposd"]["decode_osd2"] = @benchmarkable decode!(d, s) setup=(d=deepcopy($bposd2); s=copy($syn)) evals=1

# ── BP decode! ───────────────────────────────────────────────────────────────

SUITE["bp"] = BenchmarkGroup(["bp"])
bp_decoder = BeliefPropagationDecoder(H, per, 100)
SUITE["bp"]["decode"] = @benchmarkable decode!(d, s) setup=(d=deepcopy($bp_decoder); s=copy($syn)) evals=1

# ── BitFlip decode! ──────────────────────────────────────────────────────────

SUITE["bitflip"] = BenchmarkGroup(["bitflip"])
bf_decoder = BitFlipDecoder(H, per, 100)
SUITE["bitflip"]["decode"] = @benchmarkable decode!(d, s) setup=(d=deepcopy($bf_decoder); s=copy($syn)) evals=1

# ── BP-OTS decode! ───────────────────────────────────────────────────────────

SUITE["bpots"] = BenchmarkGroup(["bpots"])
bpots_decoder = BPOTSDecoder(BitMatrix(Matrix(H)), per, 100; T=9, C=2.0)
SUITE["bpots"]["decode"] = @benchmarkable decode!(d, s) setup=(d=deepcopy($bpots_decoder); s=copy($syn)) evals=1

# ── SSF (Small-Set Flip) decode! ─────────────────────────────────────────────
# Use a moderate-sized code for SSF benchmarks
H_ssf = LDPCDecoders.parity_check_matrix(100, 5, 4)
s_ssf, n_ssf = size(H_ssf)
err_ssf = rand(n_ssf) .< per
syn_ssf = Bool.((H_ssf * err_ssf) .% 2)

SUITE["ssf"] = BenchmarkGroup(["ssf"])
for t_val in (1, 2, 3)
    ssf_dec = SmallSetFlipDecoder(H_ssf, 50; t=t_val)
    SUITE["ssf"]["decode_t$t_val"] = @benchmarkable decode!(d, s) setup=(d=deepcopy($ssf_dec); s=copy($syn_ssf)) evals=1
end

# SSF batchdecode!
batch_samples = 10
errs_batch = rand(n_ssf, batch_samples) .< per
syns_batch = Bool.((H_ssf * errs_batch) .% 2)
ssf_batch_dec = SmallSetFlipDecoder(H_ssf, 50; t=2)
errs_out = zeros(Int, n_ssf, batch_samples)
conv_out = Vector{Bool}(undef, batch_samples)
SUITE["ssf"]["batchdecode"] = @benchmarkable batchdecode!(d, s, e, c) setup=(d=deepcopy($ssf_batch_dec); s=copy($syns_batch); e=copy($errs_out); c=copy($conv_out)) evals=1

# ── CSS SSF decode! ──────────────────────────────────────────────────────────
SUITE["css_ssf"] = BenchmarkGroup(["css_ssf"])
css_dec = CSSSmallSetFlipDecoder(H_ssf, H_ssf, 50; t=2)

err_x = rand(n_ssf) .< per
err_z = rand(n_ssf) .< per
syn_x = Bool.((H_ssf * err_z) .% 2)
syn_z = Bool.((H_ssf * err_x) .% 2)
syn_comb = vcat(syn_x, syn_z)

SUITE["css_ssf"]["decode_2arg"] = @benchmarkable decode!(d, sx, sz) setup=(d=deepcopy($css_dec); sx=copy($syn_x); sz=copy($syn_z)) evals=1
SUITE["css_ssf"]["decode_combined"] = @benchmarkable decode!(d, sc) setup=(d=deepcopy($css_dec); sc=copy($syn_comb)) evals=1

# CSS SSF batchdecode!
syns_comb_batch = vcat(
    Bool.((H_ssf * (rand(n_ssf, batch_samples) .< per)) .% 2),
    Bool.((H_ssf * (rand(n_ssf, batch_samples) .< per)) .% 2)
)
errs_css_out = zeros(Int, 2 * n_ssf, batch_samples)
conv_css_out = Vector{Bool}(undef, batch_samples)
SUITE["css_ssf"]["batchdecode"] = @benchmarkable batchdecode!(d, s, e, c) setup=(d=deepcopy($css_dec); s=copy($syns_comb_batch); e=copy($errs_css_out); c=copy($conv_css_out)) evals=1
