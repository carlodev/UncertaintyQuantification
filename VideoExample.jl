"""
This script is testing LegendreMultivariate testing a simple ODE, testing different parameters sensibilty: r,K

\dfrac{d X}{dt} = rX(1-\dfrac{X}{K})
"""

include("LegendreMultivariate.jl")


x0 = 70 #Initial value
dt = 0.1 #time step
timev = collect(0:dt:10)
r = 0.07
σr = 0.01

K = 210

deg = 2 #Legendre degree polynomials


#Solving the ODE
X = zeros(length(timev))
X[1] = x0

function solve_eq(x, r, K)
    r * x * (1 - x / K)
end

for (i, t) in enumerate(timev[1:end-1])
    X[i+1] = X[i] + solve_eq(X[i], r, K)
end

plot(timev, X, label="X(t)")

function solve_ode(X, r, K)
    for (i, t) in enumerate(timev[1:end-1])
        X[i+1] = X[i] + solve_eq(X[i], r, K)
    end
    return X
end

## Distributions
Ns = 1000 #sampling 

r_distribution = truncated(Normal(r, σr), r - 3 * σr, r + 3 * σr)
K_distribution = Uniform(150, 270)

kdis = create_sampling_from_distribution(K_distribution, Ns)
rdis = create_sampling_from_distribution(r_distribution, Ns)


jdis = joint_distributions(kdis,rdis)

## Compute PCE

function solve_ode(r, K)
    X = zeros(length(timev))
    X[1] = x0
    for (i, t) in enumerate(timev[1:end-1])
        X[i+1] = X[i] + solve_eq(X[i], r, K)
    end
    return X
end


## Compute K PCE

Ak = get_Legendre_matrix(kdis,deg)
bk = map(ks -> solve_ode(r, ks)[end], kdis.samples)
yk = Ak \ bk
σPCE_k = compute_σ(yk, 2)



## Compute r PCE
Ar =  get_Legendre_matrix(rdis,deg)
br = map(rs -> solve_ode(rs, K)[end], rdis.samples)
yr = Ar \ br
σPCE_r = compute_σ(yr, 2)

## Compute Joint PCE

Aj =  get_Legendre_matrix(jdis,deg)
bj = map(js -> solve_ode(js[2],js[1])[end], jdis.samples)
Aj\bj

bj_v = map(i-> map(js -> solve_ode(js[2],js[1])[i], jdis.samples),1:1:101)

yj = map(bji-> Aj\bji, bj_v)
σj = map(yji->compute_σ(yji,deg,jdis),yj)

plot(timev,X, grid=false,ribbon=σj,fillalpha=.5)

# Φ2 = map(degi -> hcubature(x->multivariate_legendre(degi, x)^2 / a_norm, -ref_domain, ref_domain)[1] ,deg_vector)
Φ2 = map(degi ->map(x->multivariate_legendre(degi, x)^2, jdis.normalized_samples) ,1:2)
