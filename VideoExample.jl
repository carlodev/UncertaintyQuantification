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
σj = map(yi->compute_σ(yi, deg),yj)

plot(timev,X, grid=false,ribbon=σj,fillalpha=.5)




function compute_σ(y, deg)
    pvector = 1:deg
    xl = collect(-1:0.0001:1) #quadrature interval
    Φ2 = map(degi -> trapz(xl, Pl.(xl, degi) .^ 2) ./ 2, pvector)
    yi = y[2:end]
    return sum(Φ2 .* yi)
end


#LHP
using LatinHypercubeSampling, QuadGK, Cuba, HCubature
plan, _ = LHCoptim(100,2,1000)
scaled_plan = scaleLHC(plan,[(-5.0,5.0),(-5.0,5.0)])
scatter(scaled_plan[:,1], scaled_plan[:,2])

i,e = hcubature(x -> cos(x[1]) * sin(x[2]), [1.0, 1.1], [2.0, 3.0])

f(x) = (Pl(x[1], 0) * Pl(x[2], 0))^2 /4
hcubature(f, [-1.0, -1.0], [1.0, 1.0])

hcubature(x-> cos(x[1]), [0.0:0.1: 1.0])
