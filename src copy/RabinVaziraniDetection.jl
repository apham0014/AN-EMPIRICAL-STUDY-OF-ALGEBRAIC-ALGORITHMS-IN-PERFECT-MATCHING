# Perfect matching detection according to Rabin & Vazirani's use of finite field arithmetic.

using Nemo
using LinearAlgebra

# This variant is our most optimised, using matrix multiplication w/ random matrix.
function DetectPerfectMatching(adjacency_matrix::Matrix{Bool})
    # Check if graph has an even number of vertices
    n = size(adjacency_matrix, 1)
    if isodd(n)
        return false
    end

    # Construct finite field (the size will be a prime that is >= n^3, where n is number of vertices in graph)
    field_size = n ^ 3
    while !is_prime(field_size)
        field_size += 1
    end
    F = GF(field_size)

    # Construct tutte matrix
    tutte_matrix = construct_tutte_matrix(adjacency_matrix, F, field_size)

    # Verify existence of perfect matching 
    determinant = det(tutte_matrix)
    if determinant == F(0)
        return false 
    else
        return true
    end
end

function construct_tutte_matrix(adjacency_matrix::Matrix{Bool}, F::FqField, field_size::Int)
    n = size(adjacency_matrix, 1)  

    # Create random integer matrix in the upper triangle
    random_ints = rand(1:(field_size-1), (n, n))
    upper_triangle = triu(random_ints, 1)
    
    # Apply mask before converting to Nemo matrix
    masked_ints = upper_triangle .* adjacency_matrix
    finite_field_matrix = matrix(F, masked_ints)

    # Make matrix skew symmetric
    tutte_matrix = finite_field_matrix - transpose(finite_field_matrix)    
    
    return tutte_matrix
end

# This variant constructs the Tutte matrix by manually creating every entry (as opposed to using matrix multiplication)
function UnoptimisedDetectPerfectMatching(adjacency_matrix::Matrix{Bool})
    # Check if graph has an even number of vertices
    n = size(adjacency_matrix, 1)
    if isodd(n)
        return false
    end

    # Construct finite field (the size will be a prime that is >= n^3, where n is number of vertices in graph)
    field_size = n ^ 3
    while !is_prime(field_size)
        field_size += 1
    end
    F = GF(field_size)

    # Construct tutte matrix
    tutte_matrix = zero_matrix(F, n, n)
    @inbounds for i in 1:n
        for j in i+1:n  
            if adjacency_matrix[i, j]
                val = rand(F)
                while val == F(0)
                    val = rand(F)
                end
                tutte_matrix[i, j] = val
                tutte_matrix[j, i] = -val
            end
        end
    end

    # Verify existence of perfect matching 
    determinant = det(tutte_matrix)
    if determinant == F(0)
        return false 
    else
        return true
    end
end