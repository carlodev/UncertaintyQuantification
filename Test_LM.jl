using Revise, Plots
include("LegendreMultivariate.jl")
deg = 2

pm = get_polynomial_combinations(deg,2) 
pol_degrees, mult_factors = find_pol_degrees(deg, pm)

X = [0.1,0.5]

multivariate_legendre(0, X::Vector; pol_degrees=pol_degrees, mult_factors=mult_factors)

#Test Video
