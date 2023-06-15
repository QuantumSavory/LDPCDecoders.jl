function base2(dec, length=0)
  binary = bitstring(UInt(dec))
  binary = split(binary, "")
  binary = [parse(Int, x) for x in binary]
  n = size(binary, 1)
  return binary[(n - length + 1):n]
end 
