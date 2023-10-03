using PolyChaos, Plots
using SparseArrays
include("UncertaintyUtils.jl")

#Example 6.1 taken from Global sensitivity analysis using polynomial chaos expansions, Bruno Sudret, https://sci-hub.hkvisa.net/10.1016/j.ress.2007.04.002

nx = 3
p = 4
Ns = 44
np = 2
Nrec = 400

Ns_theor = get_samples_number(np4, nx,p) 

Ns = Int64(Ns_theor)

model_fun(x) = prod(map(xi->3*xi^2+1,x))

uniform01 = Uniform01OrthoPoly(p; Nrec = Nrec, addQuadrature=true)
pce_uniform = convert2affinePCE(0,1,uniform01;kind="lbub")

mop = MultiOrthoPoly([uniform01, uniform01,uniform01], p)
opqp = [uniform01, uniform01,uniform01]


J_sample = PolyChaos.sampleMeasure(Ns,mop)
Ap1 = evaluate(J_sample,mop)'

bp = Float64[]


for r in 1:1:size(J_sample)[1]
    
    x_samples = evaluatePCE(pce_uniform, J_sample[r,:], uniform01)


    val = (1/2^nx) .* model_fun(x_samples)
    push!(bp,val)
end 


yp = Ap1\bp


pce_analysis(mop,yp)





## MC
mm = Float64[]
for i = 1:1:1e5
    valtmp = (1/2^nx) .*model_fun(rand(nx))
    push!(mm, valtmp)
end
mean(mm)
std(mm)




