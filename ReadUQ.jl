using Distributions, StatsPlots, Trapz
include("ReadUQ_Utils.jl")

filename = "Results/Sampling_DU89.xlsx"

res_df = read_UQ_results(filename; chord =0.2)

XLSX.writetable("DU89_results.xlsx", res_df)


x = res_df.CL ./res_df.CD
x = sort(res_df.Reat)

d = fit(Normal, x)
xn = collect(LinRange(minimum(x),maximum(x),2000))
pd = pdf.(d, xn)

plot(ylabel="PDF", xlabel="E")
plot!(xn,pd, label=false)
