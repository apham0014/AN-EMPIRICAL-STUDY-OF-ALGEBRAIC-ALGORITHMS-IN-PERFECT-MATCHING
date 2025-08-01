# Nick Harvey's perfect matching algorithm using Julia's LinearAlgebra numerical arithmetic.
# This file also includes an untested (and most-likely flawed) version of Harvey's algorithm adapted for maximum matching.

using LinearAlgebra

function FindPerfectMatching(adjacency_matrix::Matrix{Bool})
    n_original = size(adjacency_matrix, 1)
    if isodd(n_original)
        return nothing
    end

    n = nextpow(2, n_original)

    # Construct "finite field"
    field_size = n ^ 3

    # Construct tutte matrix
    tutte_matrix = zeros(Float64, n, n)
    @inbounds for i in 1:n_original
        for j in i+1:n_original  
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

    # Add edges for dummy vertices as cycle instead of complete clique
    if n_original != n
        @inbounds for i in (n_original + 1):n
            j = i+1

            if i == n
                j = n_original + 1
            end

            val = rand(1:field_size)
            while val == 0
                val = rand(1:field_size)
            end
            tutte_matrix[i,j] = val
            tutte_matrix[j,i] = -val
        end        
    end

    # Verify existence of perfect matching
    determinant = det(tutte_matrix)
    if determinant == 0
        return nothing 
    end

    # A perfect matching exists, we can therefore begin its construction
    inverse_matrix = inv(tutte_matrix)
    set_of_vertices = collect(1:n)

    # Initialize global variables
    global T = tutte_matrix
    global N = inverse_matrix

    DeleteEdgesWithin(set_of_vertices)

    # Extract matching from remaining non-zero entries in tutte
    output = Set{Tuple{Int, Int}}()
    visited = falses(n)
    for i in 1:n
        if visited[i]
            continue
        end
        for j in (i+1):n
            if T[i, j] != 0
                
                # Filter out matches with dummy edges
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

    # Base case
    n = length(S)
    if n == 1
        return
    end

    m = n ÷ 2

    # Partition S into 2 sets of the same cardinality
    S1 = S[1:m]
    S2 = S[(m+1):end]

    for (i, Si) in enumerate([S1, S2])

        Told = T[Si, Si]
        Nold = N[Si, Si]

        DeleteEdgesWithin(Si)

        # Perform rank-2 update to inverse matrix
        Delta = T[Si, Si] - Told
        if iszero(Delta)
            continue
        end

        N[Si, Si] = Nold

        # Create identity matrix of appropriate size
        T_eltype = eltype(T)
        Iblock = Matrix{T_eltype}(I, m, m)

        test_matrix = Iblock + Delta * Nold
    

        if rank(test_matrix) != size(test_matrix,1)
            # Matrix is singular - edges in Si are essential
            T[Si, Si] = Told
        else
            inv_test = inv(test_matrix)
            correction = N[S, Si] * inv_test * Delta * N[Si, S]
            N[S, S] = N[S, S] - correction

            # X = test_matrix \ Delta
            # correction = N[S, Si] * X * N[Si, S]
            # N[S, S] = N[S, S] - correction
            # N[S, S] .-= correction

        end
    end

    DeleteEdgesCrossing(S1, S2)
end

function DeleteEdgesCrossing(R::Vector{Int}, S::Vector{Int})
    global T, N
    n = length(R)

    if n == 1
        r = R[1]
        s = S[1]

        # Check if edge should be deleted
        if T[r, s] != 0 && N[r, s] != -inv(T[r, s])
            T[r, s] = 0
            T[s, r] = 0
        end
        return
    end

    m = n ÷ 2
    RS = vcat(R, S)

    # Perform paritioning
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
        eltype_T = eltype(T)
        Iblock = Matrix{eltype_T}(I, mblock, mblock)

        test_matrix = Iblock + Delta * Nold

        if rank(test_matrix) != size(test_matrix, 1)    
            # Matrix is singular - crossing edges between Ri and Sj are essential
            T[RiSj, RiSj] = Told
        else
            # Matrix is invertible - proceed with update
            inv_test = inv(test_matrix)
            correction = N[RS, RiSj] * inv_test * Delta * N[RiSj, RS]
            N[RS, RS] = N[RS, RS] - correction

            # X = test_matrix \ Delta
            # correction = N[RS, RiSj] * X * N[RiSj, RS]
            # N[RS, RS] .-= correction
        end
    end
end



# This maximum matching variant performs the same logic as the perfect matching, except on the maximum rank submatrix (which is guaranteed to have a perfect matching).
function FindMaxMatching(adjacency_matrix::Matrix{Bool})

    # Construct "finite field"
    n = size(adjacency_matrix,1)
    field_size = 0:(n ^ 3)

    # Construct tutte matrix from adjacency matrix
    tutte_matrix = zeros(n, n)
    for i in 1:n
        for j in i+1:n  
            if adjacency_matrix[i, j]
                val = rand(field_size)
                while val == 0
                    val = rand(field_size)
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

# Mucha & Sankowski's simplified algorithm for computing maximum rank submatrix 
function gaussian_elimination(A::AbstractMatrix)
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
                if A[i, j] != 0
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
