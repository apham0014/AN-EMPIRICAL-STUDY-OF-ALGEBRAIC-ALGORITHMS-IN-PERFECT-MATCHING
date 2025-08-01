# Nick Harvey's perfect matching algorithm using Nemo.jl's finite-field arithmetic.
# This file also includes an untested (and most-likely flawed) version of Harvey's algorithm adapted for maximum matching.

using Nemo

function FindPerfectMatching(adjacency_matrix::Matrix{Bool})
    # Check if graph has an even number of vertices
    n_original = size(adjacency_matrix, 1)
    if isodd(n_original)
        return nothing
    end

    # We will pad our graph later (w/ a disconnected clique of vertices) to ensure our # of vertices is a power of 2
    n = nextpow(2, n_original)

    # Construct finite field (the size will be a prime that is >= n^3, where n is number of vertices in graph)
    field_size = n ^ 3
    while !is_prime(field_size)
        field_size += 1
    end
    F = GF(field_size)

    # Construct tutte matrix
    tutte_matrix = zero_matrix(F, n, n)
    @inbounds for i in 1:n_original
        for j in i+1:n_original  
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
    
    # Add edges for dummy vertices (isolated cycle instead of complete clique)
    if n_original != n
        @inbounds for i in (n_original + 1):n
            j = i + 1

            if i == n
                j = n_original + 1
            end

            val = rand(F)
            while val == F(0)
                val = rand(F)
            end
            tutte_matrix[i,j] = val
            tutte_matrix[j,i] = -val
        end
    end

    # Verify existence of perfect matching 
    determinant = det(tutte_matrix)
    if determinant == F(0)
        return nothing 
    end

    # A perfect matching exists, we can therefore begin its construction
    inverse_matrix = inv(tutte_matrix)
    set_of_vertices = collect(1:n)
    
    # Initialize global variables
    global T = tutte_matrix
    global N = inverse_matrix
    
    DeleteEdgesWithin(set_of_vertices)

    # Return set of remaining edges
    output = Set{Tuple{Int, Int}}()
    visited = falses(n)
    for i in 1:n
        if visited[i]
            continue
        end
        for j in (i+1):n
            if T[i, j] != F(0)  

                # Filter out dummy clique matches
                if !(i > n_original && j > n_original)
                    push!(output, (i, j))
                end

                visited[i] = true
                visited[j] = true
                break
            end
        end
    end
    
    return output
end

function DeleteEdgesWithin(S::Vector{Int})
    global T, N

    # Base Case
    n = length(S)
    if n == 1
        return 
    end
    
    m = n ÷ 2
    R = base_ring(T)  
    
    # Partition into two sets
    S1 = S[1:m]
    S2 = S[(m+1):end]
    
    # Process each partition
    for (i, Si) in enumerate([S1, S2])

        Told = T[Si, Si]
        Nold = N[Si, Si]
        
        DeleteEdgesWithin(Si)

        # Perform update to inverse matrix following equation (2.3)
        Delta = T[Si, Si] - Told
        if iszero(Delta)
            continue
        end

        N[Si, Si] = Nold  
        
        Iblock = matrix_space(R, length(Si), length(Si))(1)
        test_matrix = Iblock + Delta * Nold
      
        try
            inv_test = inv(test_matrix)
            correction = N[S, Si] * inv_test * Delta * N[Si, S]
            N[S, S] = N[S, S] - correction

        catch            
            # Matrix is singular - edges in Si are essential for perfect matching
            T[Si, Si] = Told

        end
    end

    # Delete edges crossing between the two halves
    DeleteEdgesCrossing(S1, S2)
end

function DeleteEdgesCrossing(R::Vector{Int}, S::Vector{Int})
    global T, N
    F = base_ring(T)
    n = length(R)
    
    if n == 1
        r = R[1]
        s = S[1]
        
        # Check if edge should be deleted using the correct criterion
        # if T[r, s] != F(0) && N[r, s] != -inv(T[r, s])

        if T[r,s] != F(0)
            # if N[r,s] != -inv(T[r,s])
            # if N[r,s] * T[r,s] != -1
            if N[r,s] * T[r,s] != -one(F)
                T[r, s] = F(0)
                T[s, r] = F(0)
            end
        end
        return
    end

    m = n ÷ 2
    RS = vcat(R, S)
    
    # Perform partitioning
    R1 = R[1:m]
    R2 = R[(m+1):end]
    S1 = S[1:m]
    S2 = S[(m+1):end]

    for Ri in [R1, R2], Sj in [S1, S2]
        RiSj = vcat(Ri, Sj)
        Told = T[RiSj, RiSj]
        Nold = N[RiSj, RiSj]
        
        DeleteEdgesCrossing(Ri, Sj)
        
        Delta = T[RiSj, RiSj] - Told
        if iszero(Delta)
            continue
        end
        
        N[RiSj, RiSj] = Nold 
        
        mblock = length(RiSj)
        Iblock = matrix_space(F, mblock, mblock)(1)  
        test_matrix = Iblock + Delta * Nold
        
        try
            inv_test = inv(test_matrix)
            correction = N[RS, RiSj] * inv_test * Delta * N[RiSj, RS]
            N[RS, RS] = N[RS, RS] - correction
            
        catch
            # Matrix is singular - crossing edges between Ri and Sj are essential
            T[RiSj, RiSj] = Told
        end
    end
end



# This maximum matching variant performs the same logic as the perfect matching, except on the maximum rank submatrix (which is guaranteed to have a perfect matching).
function FindMaxMatching(adjacency_matrix::Matrix{Bool})

    # Construct finite field
    number_of_vertices = size(adjacency_matrix,1)
    field_size = number_of_vertices ^ 3
    while !is_prime(field_size)
        field_size += 1
    end

    F = GF(field_size)

    # Construct tutte matrix from adjacency matrix
    tutte_matrix = zero_matrix(F, number_of_vertices, number_of_vertices)
    for i in 1:number_of_vertices
        for j in i+1:number_of_vertices  
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

    # Find max rank submatrix (Ridx is ordered)
    row_pos, col_pos = gaussian_elimination(tutte_matrix)
    max_rank_submatrix = adjacency_matrix[row_pos, row_pos]
    
    submatrix_matching = FindPerfectMatching(max_rank_submatrix)

    if isnothing(submatrix_matching)
        return nothing
    end

    # Convert vertices in submatrix matching to their original values
    original_edges = Tuple{Int, Int}[]
    
    for (sub_i, sub_j) in submatrix_matching
        orig_i = row_pos[sub_i]
        orig_j = row_pos[sub_j]
        push!(original_edges, (orig_i, orig_j))
    end
    
    return original_edges
end

# Mucha & Sankowski's simplified algorithm for gaussian elimination.
function gaussian_elimination(A::FqMatrix)
    m, n = size(A)
    
    # Track which rows and columns have been used
    used_rows = falses(m)
    used_cols = falses(n)
    
    row_pos = Int[]
    col_pos = Int[]
    
    # Main elimination loop
    for step in 1:min(m, n)
        pivot_found = false
        
        # Look for a non-zero pivot in unused rows and columns
        for i in 1:m
            if used_rows[i]
                continue
            end
            
            for j in 1:n
                if used_cols[j]
                    continue
                end
                
                # Check if A[i,j] is non-zero
                if !iszero(A[i, j])
                    # Pivot found
                    push!(row_pos, i)
                    push!(col_pos, j)
                    
                    used_rows[i] = true
                    used_cols[j] = true
                    
                    pivot_found = true
                    break
                end
            end
            
            if pivot_found
                break
            end
        end
        
        if !pivot_found
            break
        end
    end
    
    return row_pos, col_pos
end
