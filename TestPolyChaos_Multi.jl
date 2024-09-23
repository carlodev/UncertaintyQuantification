using PolyChaos,Statistics,Plots
using SparseArrays
using Trapz,HCubature


x0 = 70 #Initial value
dt = 0.1 #time step
timev = collect(0:dt:10)
μr,σr  = 0.07, 0.01
μK, ΔK = 210, 60

deg_ = 2 # degree polynomials
Ns = 10000

Nrec = 100 #recursion coefficients

uniform01 = Uniform01OrthoPoly(deg_; Nrec = Nrec, addQuadrature=true)
pce_uniform = convert2affinePCE(μK-ΔK,μK+ΔK,uniform01)

mop = MultiOrthoPoly([uniform01, uniform01], deg_)



showbasis(uniform01)


J_sample = PolyChaos.sampleMeasure(Ns,mop)
Ap1 = evaluate(J_sample,mop)'

a = 0.123
b = 0.798
Jt = [a, b]
collect(evaluate(Jt,mop)')



xl = 0.0:0.00001:1.0
f0(x) = 1
f1(x) = x - 0.5
f2(x) =x^2 - 1.0x + 1/6
g0(x) = 1
g1(x) = x - 0.5
g2(x) =x^2 - 1.0x + 1/6

#This are the shape functions [OK]
f0(a) * g0(b)
f1(a)
g1(b)
f2(a)
f1(a)*g1(b)
g2(b)

ϕ0(x) = f0(x[1]) * g0(x[2])
ϕ1(x) = f1(x[1])
ϕ2(x) = g1(x[1])
ϕ3(x) = f2(x[1])
ϕ4(x) = f1(x[1])*g1(x[2])
ϕ5(x) = g2(x[2])

#Get the Tensor
_,t2= findnz(Tensor(2, mop).T)
t2
ref_domain_start = [0.0,0.0]
ref_domain_end = [1.0,1.0]

t2= Tensor(2, mop)
t2.T
t2.get([4,4]) ≈ findnz(t2.T)[2][5]

t2.op.ind

hcubature(x-> ϕ0(x).^2, ref_domain_start, ref_domain_end)[1]≈ t2[1]
hcubature(x-> ϕ1(x).^2, ref_domain_start, ref_domain_end)[1]≈ t2[2]
hcubature(x-> ϕ2(x).^2, ref_domain_start, ref_domain_end)[1]≈ t2[3]
hcubature(x-> ϕ3(x).^2, ref_domain_start, ref_domain_end)[1]≈ t2[4]
hcubature(x-> ϕ4(x).^2, ref_domain_start, ref_domain_end)[1]≈ t2[5]
hcubature(x-> ϕ5(x).^2, ref_domain_start, ref_domain_end)[1]≈ t2[6]


