using PolyChaos, Plots
using SparseArrays
using DataFrames
include("UncertaintyUtils.jl")
#Example 6.3 taken from Global sensitivity analysis using polynomial chaos expansions, Bruno Sudret, https://sci-hub.hkvisa.net/10.1016/j.ress.2007.04.002

p = 3
nx = 8
Nrec = 100
np = 10
Ns_theor = np * ((factorial(nx+p))/(factorial(nx)*factorial(p)))
Ns = Int64(ceil(Ns_theor))


sob_fun0(a,x) = (abs(4*x-2)+a)/(1+a)
sob_fun(a,x) = prod(map(sob_fun0, a, x))


a = [1,2,5,10,20,50,100,500]

uniform01 = Uniform01OrthoPoly(p; Nrec = Nrec, addQuadrature=true)
pce_uniform = convert2affinePCE(0,1,uniform01;kind="lbub")

opqp = [uniform01, uniform01,uniform01,uniform01, uniform01,uniform01,uniform01, uniform01]
mop = MultiOrthoPoly(opqp, p)

J_sample = PolyChaos.sampleMeasure(Ns,mop)
Ap1 = evaluate(J_sample,mop)'

bp = Float64[]

Samples_val = similar(J_sample)


for r in 1:1:size(J_sample)[1]
    
    x_samples = evaluatePCE(pce_uniform, J_sample[r,:], uniform01)
    Samples_val[r,:] = x_samples
    val = sob_fun(a,x_samples)
    push!(bp,val)
end 


yp = Ap1\bp

pce_analysis(mop,yp)

Statistics.cor(Samples_val[:,1],bp)
df = DataFrame(param = Samples_val[:,1], eval = bp)
df_sort = sort(df, :eval)
plot(1:1:size(df_sort)[1], df_sort.param)
plot!(1:1:size(df_sort)[1], df_sort.eval, markercolor = :orange)

scatter(1:1:size(df_sort)[1], df_sort.param, markercolor = :orange)

## MC Simulation
mm = Float64[]
for i = 1:1:1000000
    push!(mm,sob_fun(a,rand(8)))
end
mean(mm)
std(mm)