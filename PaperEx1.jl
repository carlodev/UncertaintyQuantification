using PolyChaos, Plots
using SparseArrays
include("UncertaintyUtils.jl")

#Example 6.1 taken from Global sensitivity analysis using polynomial chaos expansions, Bruno Sudret, https://sci-hub.hkvisa.net/10.1016/j.ress.2007.04.002

nx = 3 #number of uncertainty
p = 3 #polynomial degree
Ns = 29  #Number of samples
np = 2  #oversampling ratio
Nrec = 100

Ns_theor = get_samples_number(np, nx,p) 

Ns = Int64(Ns_theor)

model_fun(x) =  (1/2^nx) .* prod(map(xi->3*xi^2+1,x))

sobol_analytical(s) = 5^(-s) / ((6/5)^nx -1 )

uniform01 = Uniform01OrthoPoly(p; Nrec = Nrec, addQuadrature=true)
pce_uniform = convert2affinePCE(0,1,uniform01;kind="lbub")

mop = MultiOrthoPoly([uniform01, uniform01,uniform01], p)
opqp = [uniform01, uniform01,uniform01]

J_sample = PolyChaos.sampleMeasure(Ns,mop)

Ap1 = evaluate(J_sample,mop)'

bp = Float64[]


for r in 1:1:size(J_sample)[1]
    
    x_samples = evaluatePCE(pce_uniform, J_sample[r,:], uniform01)


    val = model_fun(x_samples)
    push!(bp,val)
end 


yp = Ap1\bp


exp_value, dev_std = pce_analysis(mop,yp)






## MC call
mm = Float64[]
for i = 1:1:1e5
    valtmp = model_fun(rand(nx))
    push!(mm, valtmp)
end

@assert isapprox(exp_value, mean(mm); rtol=0.1)
@assert isapprox(dev_std, std(mm); rtol=0.1)




