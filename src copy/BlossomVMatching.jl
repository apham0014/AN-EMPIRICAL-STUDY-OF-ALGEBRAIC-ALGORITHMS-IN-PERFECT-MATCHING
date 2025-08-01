# This file uses mlewe's Blossom V wrapper to call Kolmogorov's software.
using BlossomV

# Graph prompted to Blossom V has superfluous edges to satisfy software constraints.
function FindPerfectMatching(adjacency_matrix::Matrix{Bool})
    
    n = size(adjacency_matrix, 1)
    if isodd(n)
        return nothing
    end

    matching = Matching(n)

    # edges of original graph are added to BlossomV software as 0-cost.
    for i in 1:n
        for j in (i+1):n
            if adjacency_matrix[i, j]
                # Convert to 0-based indexing for BlossomV
                add_edge(matching, i-1, j-1, 0)
            end
        end
    end

    # add "cycle" edges to guarantee perfect matching
    for i in 1:n
        j = i + 1

        if i == n
            j = 1
        end

        if adjacency_matrix[i, j]
            continue
        else
            add_edge(matching, i-1, j-1, 1)
        end
    end

    BlossomV.solve(matching)

    # check if solution has fake/1-cost edges
    for i = 0:n-1
        paired_node = get_match(matching,i)
        if adjacency_matrix[i+1, paired_node+1]
            continue
        else
            return nothing
        end
    end

    # Convert BlossomV matrix output to set of edges.
    output = Set{Tuple{Int, Int}}()
    matches_matrix = get_all_matches(matching, n)

    length = size(matches_matrix, 2)
    for k in 1:length
        # add 1 to convert back to Julia indexing
        i = matches_matrix[1, k] + 1  
        j = matches_matrix[2, k] + 1 
        push!(output, (i, j))
    end
    
    return output
end

# Graph prompted to Blossom V has superfluous edges to satisfy software constraints.
function DetectPerfectMatching(adjacency_matrix::Matrix{Bool})
    
    n = size(adjacency_matrix, 1)
    if isodd(n)
        return false
    end

    matching = Matching(n)

    # edges of original graph are added to BlossomV software as 0-cost.
    for i in 1:n
        for j in (i+1):n
            if adjacency_matrix[i, j]
                # Convert to 0-based indexing for BlossomV
                add_edge(matching, i-1, j-1, 0)
            end
        end
    end

    # add "cycle" edges to guarantee perfect matching
    for i in 1:n
        j = i + 1

        if i == n
            j = 1
        end

        if adjacency_matrix[i, j]
            continue
        else
            add_edge(matching, i-1, j-1, 1)
        end
    end

    BlossomV.solve(matching)

    # check if solution has fake/1-cost edges
    for i = 0:n-1
        paired_node = get_match(matching,i)
        if adjacency_matrix[i+1, paired_node+1]
            continue
        else
            return false
        end
    end

    return true
end


