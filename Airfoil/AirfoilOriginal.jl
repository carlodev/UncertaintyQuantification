using PolyChaos, Plots
using SparseArrays
include("../UncertaintyUtils.jl")
include("../UQ_Airfoil_Utils.jl")

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
"""

nx = 7 #number of variables
p = 2 #polynomial degree PCE
np = 2 #oversampling ratio

Nrec = 100

Ns_theor = Int64(get_samples_number(np, nx,p))
Ns = Ns_theor 

μ0 = 1.81206e-5
rho0 = 1.225
c = 0.2
Re = 500_000
u0 = Re*μ0 / (rho0*c)

# Mean and standard deviation
U∞_mean =  37.0
U∞_σ = 0.225

TI_mean = 1.7e-3
TI_σ = 2e-4

AoA_mean = 1
AoA_σ = 0.1

μr_max = 10
μr_min = 1

σw1_mean = 0.5
σw1_σ = 0.0667

α1_mean = 0.31
α1_σ =  0.03

βstar_mean = 0.09
βstar_mean = 0.0041

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



## Multivariate
multi_dist = [U∞_dist, TI_dist,AoA_dist,μr_dist,σw1_dist,α1_dist,βstar_dist]
multi_pce = [U∞_pce, TI_pce,AoA_pce,μr_pce,σw1_pce,α1_pce,βstar_pce]

multi_mean = [U∞_mean, TI_mean, AoA_mean,μr_min, σw1_mean,α1_mean,βstar_mean ]
multi_σ = [U∞_σ,TI_σ,AoA_σ,μr_max,σw1_σ,α1_mean,βstar_mean]

var_names = ["U_freestream", "TI", "AoA", "Viscosity_Ratio", "σw1","α1","βstar"]

mop = MultiOrthoPoly(multi_dist, p)








J_sample_eval, J_sample = get_Samples(mop, Ns, multi_pce,multi_dist, multi_mean, multi_σ)

Ap1 = evaluate(J_sample,mop)'

export_samples("Sampling.xlsx", J_sample_eval, var_names)