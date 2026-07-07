using BenchmarkTools
using LDPCDecoders

const SUITE = BenchmarkGroup()

# ── Setup ────────────────────────────────────────────────────────────────────

H = LDPCDecoders.parity_check_matrix(1000, 10, 9)
per = 0.01
err = rand(1000) .< per
syn = (H * err) .% 2

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
