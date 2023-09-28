# module LegendreMultivariate
using Plots, Statistics, PolyChaos, Distributions, LegendrePolynomials
using IterativeSolvers, Trapz, Preconditioners




abstract type PCEDistribution end

struct SamplingDistribution <: PCEDistribution
    distribution
    samples::Vector{Float64}
    normalized_samples::Vector{Float64}

end

struct JointSamplingDistribution <: PCEDistribution
    distributions::Vector{SamplingDistribution}
    samples::Vector{Tuple}
    normalized_samples::Vector{Tuple}
end

function joint_distributions(d1::SamplingDistribution,d2::SamplingDistribution)
    samples = map((i,j)-> (i,j), d1.samples,d2.samples) 
    normalized_samples = map((i,j)-> (i,j), d1.normalized_samples,d2.normalized_samples) 
    JointSamplingDistribution([d1,d2], samples, normalized_samples)
end


function create_sampling_from_distribution(a_distribution, NS::Int64)
    a = minimum(a_distribution)
    b = maximum(a_distribution)
    samples = rand(a_distribution, NS)
    normalized_samples = map(s-> normalize(s, a, b),samples)
    return SamplingDistribution(a_distribution,samples,normalized_samples)

end

function get_normalized_samples(a_distribution::PCEDistribution)
    a_distribution.normalized_samples
end

a = rand(10)
b = rand(10)


hcat(a,b)


"""
    Return an array of all the tuple degrees
"""
function get_polynomial_combinations(deg::Int64,nvar::Int64)
    if nvar == 1
        pm = collect(Iterators.product(0:1:deg))

    elseif nvar == 2
        pm = collect(Iterators.product(0:1:deg,0:1:deg))

    elseif nvar == 3
        pm = collect(Iterators.product(0:1:deg,0:1:deg,0:1:deg))

    else
        @error("Degree $deg not supported, use 1,2,3")
    end
    pm = convert_iterators2Vector(pm)
    return pm
end



"""
Convert the Iterators results to a plain Vector of Tuples
"""
function convert_iterators2Vector(pm)
    a = Tuple[]

for elem in pm
    push!(a,elem)
end

return a
end



function find_pol_degrees(deg::Int64, pm::Array)
    pol_degrees = Vector[]
    mult_factor = Vector[]
    for d in 0:1:deg
        pol_degrees_tmp = Tuple[]
        mult_factor_tmp = Int64[]

        for p in pm
            if sum(p) == d
                push!(pol_degrees_tmp, p)
                fact = compute_factor( p, d)
                push!(mult_factor_tmp, fact)

            end
        end
        push!(pol_degrees, pol_degrees_tmp)
        push!(mult_factor, mult_factor_tmp)

    end
    return pol_degrees, mult_factor
end

function compute_factor(a::Tuple, order::Int64)
    factor = 1
    for i in 0:1:order
        c = count(val -> val == i, a)
        if c > 1
            factor = c*factor
        end
    end

    return 1
end






function multivariate_legendre(deg, X::Union{Tuple,Float64}; pol_degrees=nothing, mult_factors=nothing)
    nvar = length(X)

    if pol_degrees == nothing
        pm =  get_polynomial_combinations(deg,nvar)
        pol_degrees, mult_factors = find_pol_degrees(deg, pm)
    end

    res = 0.0
    for (f,d) in zip(mult_factors[deg+1],pol_degrees[deg+1])
        val_tmp = 1.0
        for (idx,xi) in zip(d,X)
            val_tmp = val_tmp * Pl(xi,idx)        
        end
        res = res + val_tmp*f
    end
return res
end


#PCE Analysis

function compute_σ(y, deg)
    pvector = 1:deg
    xl = collect(-1:0.0001:1) #quadrature interval
    Φ2 = map(degi -> trapz(xl, Pl.(xl, degi) .^ 2) ./ 2, pvector)
    yi = y[2:end]
    return sum(Φ2 .* yi)
end

"""
From arbitary interval [a;b] -> [-1,1]
"""
function normalize(val::Real, a::Real, b::Real)
    m = (a + b) * 0.5
    d = (max(a, b) - min(a, b)) * 0.5
    return (val - m) / d
end

"""
From interval  [-1,1]  -> [a;b]
"""
function scale(val::Real, a::Real, b::Real)
    m = (a + b) * 0.5
    d = (max(a, b) - min(a, b)) * 0.5
    return val * d + m
end



function get_Legendre_matrix(a_distribution::PCEDistribution,deg)
    normalized_samples = get_normalized_samples(a_distribution)
    pvector = 0:1:deg
    A1 = map(xi -> map(i -> multivariate_legendre(i, xi), pvector), normalized_samples)
    return permutedims(hcat(A1...))
     
end

# end