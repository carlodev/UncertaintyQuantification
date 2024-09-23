using PolyChaos, Plots
using SparseArrays
include("../UncertaintyUtils.jl")
include("../UQ_Airfoil_Utils.jl")
using DataFrames, XLSX

"""
Aleatoric Variables

UQ in CFD concerning
U∞ - freestream speed
TI - turbulence intensity
AoA - angle of attatck

Epistemic Variables
μr - turbulent viscoity ratio
σw1 - kw parameter
α1 - kw parameter
β* - kw parameter

γReθ model
s1 - γReθ parameter
C1 - γReθ parameter
"""

nx = 9 #number of variables
p = 2 #polynomial degree PCE
np = 2 #oversampling ratio

Nrec = 100 #value for constructing Gaussians

#number of samples
Ns_theor = Int64(get_samples_number(np, nx,p))
Ns = Int64(Ns_theor*1.2) * 10

μ0 = 1.81206e-5
rho0 = 1.225
c = 1.0 #c=0.2 for wind tunnel simulations
Re = 500_000
u0 = Re*μ0 / (rho0*c)

# Mean and standard deviation
U∞_mean =  37.0
U∞_σ = 0.225

TI_mean = 1.7e-3
TI_σ = 2e-4

AoA_mean = 1.0 #5
AoA_σ = 0.1

μr_max = 10
μr_min = 1

σw1_mean = 0.5
σw1_σ = 0.0667

α1_mean = 0.31
α1_σ =  0.03

βstar_mean = 0.09
βstar_σ  = 0.0041

s1_max = 11
s1_min = 2

C1_max = 4
C1_min = 2


#Aleatoric Variables
U∞_dist = GaussOrthoPoly(p;Nrec = Nrec, addQuadrature=true);
U∞_pce = convert2affinePCE(U∞_mean, U∞_σ, U∞_dist)

TI_dist = GaussOrthoPoly(p;Nrec = Nrec, addQuadrature=true);
TI_pce = convert2affinePCE(TI_mean, TI_σ, TI_dist)

AoA_dist = GaussOrthoPoly(p;Nrec = Nrec, addQuadrature=true);
AoA_pce = convert2affinePCE(AoA_mean, AoA_σ, AoA_dist)


#Epistemic Variables

μr_dist = Uniform01OrthoPoly(p; Nrec = Nrec, addQuadrature=true)
μr_pce = convert2affinePCE(μr_min,μr_max,μr_dist;kind="lbub")

σw1_dist = GaussOrthoPoly(p;Nrec = Nrec, addQuadrature=true);
σw1_pce = convert2affinePCE(σw1_mean, σw1_σ, σw1_dist)

α1_dist = GaussOrthoPoly(p;Nrec = Nrec, addQuadrature=true);
α1_pce = convert2affinePCE(α1_mean,α1_σ, α1_dist)

βstar_dist = GaussOrthoPoly(p;Nrec = Nrec, addQuadrature=true);
βstar_pce = convert2affinePCE(βstar_mean, βstar_σ, βstar_dist)

s1_dist = Uniform01OrthoPoly(p; Nrec = Nrec, addQuadrature=true)
s1_pce = convert2affinePCE(s1_min,s1_max,s1_dist;kind="lbub")

C1_dist = Uniform01OrthoPoly(p; Nrec = Nrec, addQuadrature=true)
C1_pce = convert2affinePCE(C1_min,C1_max,C1_dist;kind="lbub")


## Multivariate
multi_dist = [U∞_dist, TI_dist,AoA_dist,μr_dist,σw1_dist,α1_dist,βstar_dist,s1_dist,C1_dist]
multi_pce = [U∞_pce, TI_pce,AoA_pce,μr_pce,σw1_pce,α1_pce,βstar_pce,s1_pce,C1_pce]

multi_mean = [U∞_mean, TI_mean, AoA_mean,μr_min, σw1_mean,α1_mean,βstar_mean,s1_min,C1_min ]
multi_σ = [U∞_σ,TI_σ,AoA_σ,μr_max,σw1_σ,α1_σ,βstar_σ,s1_max,C1_max]

var_names = ["U_freestream", "TI", "AoA", "Viscosity_Ratio", "σw1","α1","βstar","s1", "C1"]

mop = MultiOrthoPoly(multi_dist, p)

using JLD2
save("MOP_DU89_$(Re)_$(AoA_mean).jld2", "mop", mop)





J_sample_eval, J_sample = get_Samples(mop, Ns, multi_pce,multi_dist, multi_mean, multi_σ; idxs=[1,2,3,5,6,7], half_normal=6)



Ap1 = evaluate(J_sample,mop)'

export_samples("Sampling_DU89_$(Re)_$(AoA_mean).xlsx", J_sample_eval, var_names)

export_samples("SamplingN_DU89_$(Re)_$(AoA_mean).xlsx", J_sample, var_names)



Ap1 = evaluate(J_sample,mop)'



histogram(J_sample_eval[:,6])
histogram(J_sample[:,6])


## Gaussian mean 0.0, std 1.0; always excluding values larger than ±3 
mean(J_sample[:,1])
std(J_sample[:,1])

mean(J_sample[:,2])
std(J_sample[:,2])


mean(J_sample[:,3])
std(J_sample[:,3])


## mean half gaussian
sqrt(2/pi)
mean(J_sample[:,6])

## std half gaussian
sqrt((1-2/pi))
std(J_sample[:,6])

