function base2(dec, length=0)
  binary = bitstring(UInt(dec))
  binary = split(binary, "")
  binary = [parse(Int, x) for x in binary]
  n = size(binary, 1)
  return binary[(n - length + 1):n]
end 

function row_echelon(matrix, full=false)
    numrows, numcols = size(matrix)
    pivot_row = 1
    pivot_col = []
    
    the_matrix = matrix
    for col in range(1,numcols)
          if the_matrix[pivot_row, col] != 1
              swap_row_index = pivot_row + findmax(the_matrix[pivot_row:numrows, col])

              if the_matrix[swap_row_index, col] == 1
                the_matrix[[swap_row_index, pivot_row]] = the_matrix[[pivot_row, swap_row_index]]
                
                transform_matrix[[swap_row_index, pivot_row]] = transform_matrix[[pivot_row, swap_row_index]]
              end

          end
          if the_matrix[pivot_row, col]

            if not full
                elimination_range = [k for k in range(pivot_row + 1, num_rows)]
            else
                elimination_range = [k for k in range(1:num_rows) if k != pivot_row]
            end

            for j in elimination_range

              if the_matrix[j, col] != 0 && pivot_row != j

                  the_matrix[j] = (the_matrix[j] + the_matrix[pivot_row]) % 2

                  transform_matrix[j] = (transform_matrix[j] + transform_matrix[pivot_row]) % 2
              end
            end
            pivot_row += 1
            pivot_cols.append(col)
          end

        if pivot_row >= numrows
            break
        end
      end
    # The rank is equal to the maximum pivot index
    matrix_rank = pivot_row
    row_esch_matrix = the_matrix

    return [row_esch_matrix, matrix_rank, transform_matrix, pivot_cols]
end


