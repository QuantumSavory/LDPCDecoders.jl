# Optional test flags
JET_flag = ARGS == ["jet"]

if JET_flag
  @info "Running JET tests in their dedicated test environment."
  using Pkg
  Pkg.activate(joinpath(@__DIR__, "projects", "jet"))
  Pkg.instantiate()
else
  @info "Skipping JET tests -- pass `test_args=[\"jet\"]` to Pkg.test to enable them."
end

using TestItemRunner
using LDPCDecoders

# filter for the test
testfilter = ti -> begin
  exclude = Symbol[]

  if JET_flag
    return :jet in ti.tags
  else
    push!(exclude, :jet)
  end

  if !(VERSION >= v"1.10")
    push!(exclude, :doctests)
    push!(exclude, :aqua)
  end

  return all(!in(exclude), ti.tags)
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")
println("ENV[\"PYTHON\"] = \"$(get(ENV,"PYTHON",nothing))\"")

@run_package_tests filter=testfilter
