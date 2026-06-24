@testitem "Aqua" tags=[:aqua] begin
using Aqua
using LDPCDecoders
Aqua.test_all(LDPCDecoders)
end
