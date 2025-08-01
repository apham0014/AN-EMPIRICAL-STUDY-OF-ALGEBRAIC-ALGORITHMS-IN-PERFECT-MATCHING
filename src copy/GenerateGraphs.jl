# This file contains my code for creating sample graphs. 
# All functions uniquely output files to my personal file directories.

# My generated graphs are text files with the following format:
#   - First line of the .txt file is the number of vertices.
#   - Every succeeding line is an edge denoted by "i, j".

using Random

function generateCompleteGraph(n::Integer)
    fname = "K$(n).txt"
    dirpath = joinpath(homedir(), "Documents", "DSRI", "graphs")    
    fpath = joinpath(dirpath, fname)

    open(fpath, "w") do file
        # first line is number of vertices
        println(file, n)

        # succeeding lines are edges
        for i = 1:n
            for j = i+1:n
                if i != j
                    println(file, "$i, $j")
                end
            end
        end
    end
    println("Generated successfully.")

end

function generateCycleGraph(n::Integer)
    fname = "C$(n).txt"
    dirpath = joinpath(homedir(), "Documents", "DSRI", "graphs")    
    fpath = joinpath(dirpath, fname)

    open(fpath, "w") do file
        # first line is number of vertices
        println(file, n)

        # succeeding lines are edges
        for i = 1:n
            j = i + 1
            if i < n
                println(file, "$i, $j")
            else
                println(file, "$i, 1")
            end
        end
    end
    println("Generated successfully.")
end

# Perfect matching exists. Graph density > 10%
function generateSparseGraph(n::Int)
    if n % 2 != 0
        return   
    end
    
    fname = "Sparse$(n).txt"
    dirpath = joinpath(homedir(), "Documents", "DSRI", "graphs")    
    fpath = joinpath(dirpath, fname)

    open(fpath, "w") do file
        # first line is number of vertices
        println(file, n)


        # add edges to guarantee perfect matching
        adj = falses(n, n)
        vertices = collect(1:n)
        shuffle!(vertices)

        for i in 1:2:n-1
            v1 = vertices[i]
            v2 = vertices[i+1]
            println(file, "$v1, $v2")
            adj[v1 ,v2] = true
            adj[v2, v1] = true
        end


        # add superfluous edges 
        for v1 = 1:n
            for v2 = v1+1:n
            
                # 10 percent probability
                should_add_edge = rand(1:10)
                if should_add_edge == 1 && !adj[v1, v2]
                    println(file, "$v1, $v2")
                    adj[v1, v2] = adj[v2, v1] = true
                end
            end
        end
    end
    println("Generated successfully.")
end

# Perfect matching exists. Graph density > 70%
function generateDenseGraph(n::Int)
    if n % 2 != 0
        return   
    end
    
    fname = "Dense$(n).txt"
    dirpath = joinpath(homedir(), "Documents", "DSRI", "graphs")    
    fpath = joinpath(dirpath, fname)

    open(fpath, "w") do file
        # first line is number of vertices
        println(file, n)


        # add edges for perfect matching
        adj = falses(n, n)
        vertices = collect(1:n)
        shuffle!(vertices)

        for i in 1:2:n-1
            v1 = vertices[i]
            v2 = vertices[i+1]
            println(file, "$v1, $v2")
            adj[v1 ,v2] = true
            adj[v2, v1] = true
        end


        # add superfluous edges 
        for v1 = 1:n
            for v2 = v1+1:n
            
                # 70 percent probability
                should_add_edge = rand(1:10)
                if should_add_edge < 8 && !adj[v1, v2]
                    println(file, "$v1, $v2")
                    adj[v1, v2] = adj[v2, v1] = true
                end
            end
        end
    end
    println("Generated successfully.")
end

# No perfect matching because two vertices are isolated
function generateNonPerfectGraph(n::Int)
    if n % 2 != 0
        println("Error: n must be even")
        return   
    end
    
    fname = "NonPerfect$(n).txt"
    dirpath = joinpath(homedir(), "Documents", "DSRI", "graphs")    
    fpath = joinpath(dirpath, fname)

    open(fpath, "w") do file
        println(file, n)

        # Create adjacency matrix
        adj = falses(n, n)
        vertices = collect(1:n)
        shuffle!(vertices)
        
        # Remove isolated vertices
        non_isolated = vertices[1:end-2]  
        
        for i in 1:2:length(non_isolated)-1  
            v1 = non_isolated[i]
            v2 = non_isolated[i+1]
            println(file, "$v1, $v2")
            adj[v1, v2] = true
            adj[v2, v1] = true
        end

        # Add superfluous edges among non-isolated vertices
        for v1 in non_isolated
            for v2 in non_isolated
                if v1 < v2 
                    # 50% probability of adding edge if not already present
                    should_add_edge = rand(1:10)
                    if should_add_edge < 6 && !adj[v1, v2]
                        println(file, "$v1, $v2")
                        adj[v1, v2] = true
                        adj[v2, v1] = true
                    end
                end
            end
        end
    end
    println("Generated successfully.")
end

# No perfect matchings in star graphs.
function generateStarGraph(n::Int)
    if n % 2 != 0
        return   
    end
    
    fname = "Star$(n).txt"
    dirpath = joinpath(homedir(), "Documents", "DSRI", "graphs")    
    fpath = joinpath(dirpath, fname)

    open(fpath, "w") do file
        # first line is number of vertices
        println(file, n)

        for i in 2:n
            println(file, "1, $(i)")
        end
    end
    println("Generated successfully.")
end

# No perfect matching via tutte's condition: our graph will have one subset S of V where (odd components of G\S) > cardinality of S.
function generateNonPerfectGraphB(n::Int)
    if n % 2 != 0
        return   
    end
    
    fname = "NonPerfectB$(n).txt"
    dirpath = joinpath(homedir(), "Documents", "DSRI", "graphs")    
    fpath = joinpath(dirpath, fname)

    open(fpath, "w") do file
        # first line is number of vertices
        println(file, n)

        # our graph has 5 components in the case: V \ {1}
        component_1 = Vector{Int}()
        component_2 = Vector{Int}()
        component_3 = Vector{Int}()
        component_4 = Vector{Int}()
        component_5 = Vector{Int}()
        components = [component_1, component_2, component_3, component_4, component_5]

        for i = 2:6
            println(file, "1, $i")
            push!(components[i-1], i)
        end

        n_current = 7

        # we will add 2 new vertices at a time, to maintain odd cardinality of components
        while n_current < n

            # connect 2 new vertices together
            println(file, "$(n_current), $(n_current + 1)")

            # which component do we add to?
            decision_1 = rand(1:5)

            # connect both new vertices to component
            for i = 1:2
                
                # how many vertices should we connect to?
                decision_2 = rand(1:length(components[decision_1]))
                
                # connect to most recent vertices in selected component
                for j in 1:decision_2
                    v2 = components[decision_1][length(components[decision_1]) - j + 1]
                    println(file, "$n_current, $v2")
                end

                n_current += 1
            end
        end
    end

    println("Generated successfully.")
end


