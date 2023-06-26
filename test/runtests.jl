using LDPC
using Test

@testset "LDPC.jl" begin
    
    error_rate = 0.1
    
    parity_check_matrix = [1 0 1 0 1 0 1; 0 1 1 0 0 1 1; 0 0 0 1 1 1 1]
    actual_message = [1 1 1 1 1 1 1]
    received_message = [1 1 1 1 1 0 1]  # 2nd bit from right is an error
    decoded_message, success = bp_decode(parity_check_matrix, received_message, error_rate)
    
    actual_message = vec(actual_message)
    @test decoded_message == actual_message

    parity_check_matrix = [1 0 1 1 1 0 0; 1 1 0 1 0 1 0; 1 1 1 0 0 0 1]
    actual_message = [1 0 0 0 1 1 1]
    received_message = [1 0 0 1 1 1 1] # 4th bit from left is an error
    decoded_message, success = bp_decode(parity_check_matrix, received_message, error_rate)
    
    actual_message = vec(actual_message)
    @test decoded_message == actual_message
    
    parity_check_matrix = [1 1 1 1 0 0 0 0 0 0; 1 0 0 0 1 1 1 0 0 0; 0 1 0 0 1 0 0 1 1 0; 0 0 1 0 0 1 0 1 0 1; 0 0 0 1 0 0 1 0 1 1]
    actual_message = [1 1 0 0 1 0 0 0 0 0]
    received_message = [1 1 0 0 1 0 0 0 0 1] # last bit is an error
    decoded_message, success = bp_decode(parity_check_matrix, received_message, error_rate)
    
    actual_message = vec(actual_message)
    @test decoded_message == actual_message

end
