using SafeTestsets
using LDPCDecoders

function doset(descr)
    if length(ARGS) == 0
        return true
    end
    for a in ARGS
        if occursin(lowercase(a), lowercase(descr))
            return true
        end
    end
    return false
end

macro doset(descr)
    quote
        if doset($descr)
            @safetestset $descr begin
                include("test_"*$descr*".jl")
            end
        end
    end
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")
println("ENV[\"PYTHON\"] = \"$(get(ENV,"PYTHON",nothing))\"")



@doset "oldtests"
@doset "bp_decoder"
@doset "bf_decoder"

VERSION >= v"1.10" && @doset "doctests"
get(ENV,"JET_TEST","")=="true" && @doset "jet"
VERSION >= v"1.10" && @doset "aqua"
