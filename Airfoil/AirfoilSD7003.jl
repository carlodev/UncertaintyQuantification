using PolyChaos, Plots
using SparseArrays
using DataFrames, XLSX
include("../UncertaintyUtils.jl")
include("../UQ_Airfoil_Utils.jl")

"""
Aleatoric Variables

UQ in CFD concerning
TI - turbulence intensity

Epistemic Variables

μr - turbulent viscoity ratio
σw1 - kw parameter
α1 - kw parameter
β* - kw parameter

γReθ model
s1 - γReθ parameter
C1 - γReθ parameter
"""

nx = 7 #number of variables
p = 2 #polynomial degree PCE
np = 2 #oversampling ratio

Nrec = 100

Ns_theor = Int64(get_samples_number(np, nx,p))

Ns = Ns_theor 

Re = 60_000
u0 = 1.0
rho0 = 1.0
c = 1.0
μ0 = (rho0*c) /(Re)

# Mean and standard deviation

TI_mean = 5.0e-4
TI_σ = 2e-4


μr_max = 10
μr_min = 1

σw1_mean = 0.5
σw1_σ = 0.0667

α1_mean = 0.31
α1_σ =  0.03

βstar_mean = 0.09
βstar_σ = 0.0041

s1_max = 11
s1_min = 2

C1_max = 4
C1_min = 2


#Aleatoric Variables

TI_dist = GaussOrthoPoly(p;Nrec = Nrec, addQuadrature=true);
TI_pce = convert2affinePCE(TI_mean, TI_σ, TI_dist)


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
multi_dist = [TI_dist,μr_dist,σw1_dist,α1_dist,βstar_dist,s1_dist,C1_dist]
multi_pce = [TI_pce,μr_pce,σw1_pce,α1_pce,βstar_pce,s1_pce,C1_pce]

multi_mean = [TI_mean, μr_min, σw1_mean,α1_mean,βstar_mean,s1_min,C1_min ]
multi_σ = [TI_σ,μr_max,σw1_σ,α1_mean,βstar_mean,s1_max,C1_max]


var_names = ["TI", "Viscosity_Ratio", "σw1","α1","βstar", "s1", "C1"]

mop = MultiOrthoPoly(multi_dist, p)




J_sample_eval, J_sample = get_Samples(mop, Ns, multi_pce, multi_dist, multi_mean, multi_σ; idxs=[1,3,4,5], half_normal=4)

Ap1 = evaluate(J_sample,mop)'


export_samples("Sampling_sd7003.xlsx", J_sample_eval, var_names)
export_samples("SamplingN_sd7003.xlsx", J_sample, var_names)

