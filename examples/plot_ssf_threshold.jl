# SSF Threshold Plotting Script
#
# Reads the CSV output from ssf_threshold.jl and generates waterfall
# threshold curves using CairoMakie. Follows the same plotting pattern
# used in QuantumClifford.jl's ECC evaluation docs.
#
# Prerequisites:
#   ] add CairoMakie DelimitedFiles
#
# Usage:
#   julia examples/plot_ssf_threshold.jl                          # uses default file
#   julia examples/plot_ssf_threshold.jl my_results.csv           # custom CSV
#   julia examples/plot_ssf_threshold.jl results.csv out.png      # custom output

using CairoMakie
using DelimitedFiles
using Printf

# ============================================================================
# Read CSV
# ============================================================================

csv_file = length(ARGS) >= 1 ? ARGS[1] : "ssf_threshold_results.csv"
out_file = length(ARGS) >= 2 ? ARGS[2] : replace(csv_file, ".csv" => ".png")

if !isfile(csv_file)
    error("CSV file not found: $csv_file\nRun `julia examples/ssf_threshold.jl` first to generate it.")
end

raw = readdlm(csv_file, ','; header=true)
data = raw[1]
header = raw[2]

distances    = Int.(data[:, 1])
p_values     = Float64.(data[:, 2])
logical_err  = Float64.(data[:, 3])
decode_fail  = Float64.(data[:, 4])

unique_L = sort(unique(distances))

# ============================================================================
# Threshold crossing estimation (linear interpolation)
# ============================================================================

function find_crossings(L1, L2, distances, p_values, logical_err)
    mask1 = distances .== L1
    mask2 = distances .== L2
    p1, ler1 = p_values[mask1], logical_err[mask1]
    p2, ler2 = p_values[mask2], logical_err[mask2]

    crossings = Float64[]
    n = min(length(ler1), length(ler2))
    for j in 1:n-1
        d_j  = ler2[j]   - ler1[j]
        d_j1 = ler2[j+1] - ler1[j+1]
        if d_j * d_j1 < 0  # sign change
            frac = abs(d_j) / (abs(d_j) + abs(d_j1))
            push!(crossings, p1[j] + frac * (p1[j+1] - p1[j]))
        end
    end
    return crossings
end

# estimate threshold from all adjacent pairs
all_crossings = Float64[]
for i in 1:length(unique_L)-1
    append!(all_crossings, find_crossings(unique_L[i], unique_L[i+1],
                                           distances, p_values, logical_err))
end
p_th = isempty(all_crossings) ? nothing : sum(all_crossings) / length(all_crossings)

# ============================================================================
# Plot: Logical Error Rate waterfall
# ============================================================================

colors = Makie.wong_colors()

fig = Figure(size = (700, 500), fontsize = 14)
ax = Axis(fig[1, 1],
    xlabel = "Physical Error Rate (p)",
    ylabel = "Logical Error Rate (P_L)",
    title  = "SSF Decoder Threshold Curves (Toric Code HGP)",
)

for (idx, L) in enumerate(unique_L)
    mask = distances .== L
    scatterlines!(ax, p_values[mask], logical_err[mask];
        color = colors[idx],
        marker = :circle,
        markersize = 8,
        linewidth = 2,
        label = "L=$L",
    )
end

# threshold line
if p_th !== nothing
    vlines!(ax, [p_th]; color = :gray, linestyle = :dash, linewidth = 1.5,
            label = @sprintf("p_th ≈ %.1f%%", p_th * 100))
end

axislegend(ax; position = :rb, framevisible = true, labelsize = 12)

save(out_file, fig, px_per_unit = 2)
println("Saved: $out_file")

# ============================================================================
# Plot: Decode Failure Rate (separate figure)
# ============================================================================

fig2 = Figure(size = (700, 500), fontsize = 14)
ax2 = Axis(fig2[1, 1],
    xlabel = "Physical Error Rate (p)",
    ylabel = "Decode Failure Rate",
    title  = "SSF Decode Failure Rate (Toric Code HGP)",
)

for (idx, L) in enumerate(unique_L)
    mask = distances .== L
    scatterlines!(ax2, p_values[mask], decode_fail[mask];
        color = colors[idx],
        marker = :diamond,
        markersize = 8,
        linewidth = 2,
        label = "L=$L",
    )
end

axislegend(ax2; position = :lt, framevisible = true, labelsize = 12)

fail_out = replace(out_file, ".png" => "_decode_failures.png")
save(fail_out, fig2, px_per_unit = 2)
println("Saved: $fail_out")
