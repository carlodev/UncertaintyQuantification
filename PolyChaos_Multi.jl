using PolyChaos,Statistics,Plots
using SparseArrays

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


muv = [μK-ΔK, μr]
sigv = [μK+ΔK,σr ]
opqv = [uniform01, gaussian01]

fun = (x,y,z)->convert2affinePCE(x, y, z)
pce_mop1 = map(fun, muv, sigv, opqv)
pce_mop2= map(i-> assign2multi(pce_mop1[i], i, mop.ind), 1:2)
pce_mop = [findnz(pce_mop2[1])[2]...,findnz(pce_mop2[2])[2]...]



t2 = Tensor(2,mop)
_, T2 = findnz(t2.T)
t3 = Tensor(3,mop)
_, T3 = findnz(t3.T)

J_sample = PolyChaos.sampleMeasure(Ns,mop)
k_samples_unifrom = evaluatePCE(pce_uniform, J_sample[:,1],uniform01)
r_samples_gaussian = evaluatePCE(pce_gaussian, J_sample[:,2],gaussian01)



map_sol = (r,k)->solve_ode(r,k)[end]

bk = map(map_sol, r_samples_gaussian,k_samples_unifrom)

Aj1 = evaluate(J_sample,mop)'

yk = Aj1\bk

yk[1]
D = sum(T2[2:end] .* yk[2:end].^2 )
σ =  sqrt(D) 

(T2[2:end] .* yk[2:end].^2) ./ D

D_K = sum(T2[2] .* yk[2].^2  + T2[2] .* yk[3].^2  )
D_r = sum(T2[3] .* yk[3].^2 )

T2[2] .* yk[3].^2
T2[3] .* yk[2].^2







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


