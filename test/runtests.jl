# Optional test flags
JET_flag = get(ENV, "JET_TEST", "") == "true"

JET_flag && @info "Running with JET tests."

using Pkg
JET_flag && Pkg.add("JET")

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
