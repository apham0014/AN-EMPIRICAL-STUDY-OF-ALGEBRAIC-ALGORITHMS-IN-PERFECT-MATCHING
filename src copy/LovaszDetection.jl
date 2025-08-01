# Perfect matching detection according to Lovasz's use of the rationals.

using LinearAlgebra


# Detection process constructs Tutte matrix with floating point values.
function FloatDetectPerfectMatching1(adjacency_matrix::Matrix{Bool})
    # Check if graph has an even number of vertices
    n = size(adjacency_matrix, 1)
    if isodd(n)
        return false
    end

    # Construct tutte matrix
    tutte_matrix = construct_tutte_matrix(adjacency_matrix)

    # Verify existence of perfect matching
    determinant = det(tutte_matrix)
    if abs(determinant) < 1e-10
        return false 
    else
        return true
    end
end

function construct_tutte_matrix(adjacency_matrix::Matrix{Bool})
    n = size(adjacency_matrix, 1)  
    
    # Generate square matrix of random values
    random_vals = rand(Float64, (n, n))

    
    # Extract upper triangle part
    upper_tri = triu(random_vals, 1)
    
    # Apply mask
    masked = upper_tri .* adjacency_matrix  

    # Make matrix skew symmetric
    tutte_matrix = masked - masked'
    
    return tutte_matrix
end

# Detection process constructs Tutte matrix with Int values.
function IntDetectPerfectMatching(adjacency_matrix::Matrix{Bool})
    # Check if graph has an even number of vertices
    n = size(adjacency_matrix, 1)
    if isodd(n)
        return false
    end

    # Construct "finite field"
    field_size = n ^ 3

    # Construct tutte matrix
    tutte_matrix = zeros(Int, n, n)
    @inbounds for i in 1:n
        for j in i+1:n  
            if adjacency_matrix[i, j]
                val = rand(1:field_size)
                while val == 0
                    val = rand(1:field_size)
                end
                tutte_matrix[i, j] = val
                tutte_matrix[j, i] = -val
            end
        end
    end

    # Verify existence of perfect matching
    determinant = det(tutte_matrix)
    if abs(determinant) < 1e-10        
        return false 
    else
        return true
    end
end