using PolyChaos, SparseArrays, Combinatorics

using Statistics, Trapz,HCubature, Plots, LinearAlgebra



function compute_Sobolov_index(v::Vector, t2::Tensor, Di::Vector)

    @assert length(findnz(t2.T)[2]) == length(Di)

    indexes = t2.op.ind #Get the matrix of Shape fun indexes
    nvar = length(t2.op.uni) #Get the number of indipendent variables
    
    index_total = collect(1:1:nvar)
    
    idx_not_empty = Vector{Int64}[] #where column of v not 0
    for vi in v
        idx_tmp = findall(x-> x > 0, indexes[:,vi])
        push!(idx_not_empty,idx_tmp)
    end
    
    idx_empty = Vector{Int64}[] #where column of not v are 0
    for vi in setdiff(index_total,v)
        idx_tmp = findall(x-> x == 0, indexes[:,vi])
        push!(idx_empty,idx_tmp)
    end
    
    
    idx_v = intersect(idx_empty...,idx_not_empty...)
    
    # println(t2.op.ind[idx_v,:])
    
    return sum(Di[idx_v])
    end
    
    
    function compute_Sobolov_index(v::Int64, t2::Tensor, Di::Vector)
        compute_Sobolov_index([v], t2::Tensor, Di::Vector)
    end


    function compute_multiple_Sobolov_index(nvar::Int64,t2::Tensor, Di::Vector )
    

        comb = Vector{Vector{Int64}}[]
        for a = 1:1:nvar
            comb_tmp = collect(combinations(1:nvar,a))
            push!(comb, comb_tmp)
        end
        
        for idx_order in comb[1:2] #First and second order indexes
            for idx_ in idx_order
                s_tmp = compute_Sobolov_index(idx_, t2, Di)
                println("S$(idx_): $(s_tmp)")
            end
        end
    end

    
function pce_analysis(mop::AbstractOrthoPoly, y::Vector)
    t2 = Tensor(2,mop) #Tensor 2
    _,T2 = findnz(t2.T)
    mop_deg = mop.deg

    println("PCE Analysis")
    println("Polynomial degree: $(mop_deg)")
    nvar = length(t2.op.uni) #Get the number of indipendent variables
    println("Number of independet variables: $(nvar)")
    println("Total number of PCE coefficients: $(length(y))")



    expected_value = y[1]
    println("Expected Value: $(expected_value)")

 
    
    D = sum(T2[2:end] .* y[2:end].^2 ) #variance
    println("Variance: $(D)")
    σ = sqrt(D)
    println("Deviation standard (σ): $( σ)")
    
    if nvar>1
        Di =  (T2.* y.^2) ./ D
        compute_multiple_Sobolov_index(nvar,t2, Di )
    end

    return expected_value, σ
end

"""
    get_samples_number(np::Int64, nx::Int64,p::Int64)

np: oversampling ratio (np = 2)
nx:: number of independent variables
p:: polynomial degree
"""
function get_samples_number(np::Int64, nx::Int64,p::Int64)
    return np * ((factorial(nx+p))/(factorial(nx)*factorial(p)))
end
