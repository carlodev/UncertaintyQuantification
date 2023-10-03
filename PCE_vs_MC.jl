using PolyChaos,Statistics,Plots
using SparseArrays
using IterativeSolvers
using LaTeXStrings
include("UncertaintyUtils.jl")


function solve_eq(x, r, K)
    r * x * (1 - x / K)
 end
 
function solve_ode(r, K)
    X = zeros(length(timev))
    X[1] = x0
    for (i, t) in enumerate(timev[1:end-1])
        X[i+1] = X[i] + solve_eq(X[i], r, K)
    end
    return X
end




x0 = 70 #Initial value
dt = 0.1 #time step
timev = collect(0:dt:10)
μr,σr  = 0.07, 0.01
μK, ΔK = 210, 60

deg_ = 2 # degree polynomials
Ns = 10000

Nrec = 40 #recursion coefficients


uniform01 = Uniform01OrthoPoly(deg_; Nrec = Nrec, addQuadrature=true)
pce_uniform = convert2affinePCE(μK-ΔK,μK+ΔK,uniform01;kind="lbub")

gaussian01 = GaussOrthoPoly(deg_;Nrec = Nrec, addQuadrature=true);
pce_gaussian = convert2affinePCE(μr,σr,gaussian01)

mop = MultiOrthoPoly([uniform01, gaussian01], deg_)





mean_pce = Float64[]
σ_pce = Float64[]

mean_mc = Float64[]
σ_mc = Float64[]


ns_vector = [collect(2:9)... , 10 .^(collect(1:1:6))...,12]

np = 2
Ns_theor = np * ((factorial(2+2))/(factorial(2)*factorial(2)))

for ns in ns_vector
    J_sample = PolyChaos.sampleMeasure(ns,mop)
    k_samples_unifrom = evaluatePCE(pce_uniform, J_sample[:,1],uniform01)
    r_samples_gaussian = evaluatePCE(pce_gaussian, J_sample[:,2],gaussian01)
    map_sol = (r,k)->solve_ode(r,k)[end]

    bk = map(map_sol, r_samples_gaussian,k_samples_unifrom)
    
    Aj1 = evaluate(J_sample,mop)'
    yk = ones(size(Aj1)[2])
    
    
    yk = Aj1\bk

    ev,σ =  pce_analysis(mop,yk)
    
    push!(mean_pce,ev)
    push!(σ_pce,σ)

    push!(mean_mc, mean(bk))
    push!(σ_mc, std(bk))


end



plt_mean = plot(title = "Deviation Standard value convergence", xlabel = "Number of samples", ylabel ="expected value" )
plot!(ns_vector[1:end-1], mean_pce[1:end-1], xaxis=:log, label = "PCE",xticks = ([0,10,100,1000,10000,100000,1e6], [latexstring("\$0\$"), 
latexstring("\$10\$"), latexstring("\$10^2\$"), latexstring("\$10^3\$"),latexstring("\$10^4\$"),latexstring("\$10^5\$"),latexstring("\$10^6\$")]))
plot!(ns_vector[1:end-1],mean_mc[1:end-1], xaxis=:log, label = "MC simulation" )
scatter!([ns_vector[end]], [mean_pce[end]], xaxis=:log, label = "PCE - np=2", markershape = :circle)
savefig(plt_mean, "MeanValue_conv.pdf")


plt_std = plot(title = "Deviation Standard value convergence", xlabel = "Number of samples", ylabel =latexstring("\$\\sigma\$") )
plot!(ns_vector[1:end-1], σ_pce[1:end-1], xaxis=:log, label = "PCE",xticks = ([0,10,100,1000,10000,100000,1e6], [latexstring("\$0\$"), 
latexstring("\$10\$"), latexstring("\$10^2\$"), latexstring("\$10^3\$"),latexstring("\$10^4\$"),latexstring("\$10^5\$"),latexstring("\$10^6\$")]))
plot!(ns_vector[1:end-1],σ_mc[1:end-1], xaxis=:log, label = "MC simulation" )
scatter!([ns_vector[end]], [σ_pce[end]], xaxis=:log, label = "PCE - np=2", markershape = :circle)
savefig(plt_std, "Std_conv.pdf")




