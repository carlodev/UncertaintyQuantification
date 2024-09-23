using Distributions, XLSX, DataFrames, Distributions, StatsPlots, CategoricalArrays,Statistics
using  IterativeSolvers, StatsBase

include("../UncertaintyUtils.jl")

include("../ReadUQ_Utils.jl")

airfoilname = "DU89"
sampling_df = DataFrame(XLSX.readtable("Airfoil/Sampling_$(airfoilname).xlsx", "Sheet1"))
sampling_dfN = DataFrame(XLSX.readtable("Airfoil/SamplingN_$(airfoilname).xlsx", "Sheet1"))

minimum(sampling_df.α1)
std(sampling_df.α1)

res_df = DataFrame(XLSX.readtable("Airfoil/Results/$(airfoilname)_results.xlsx", "Sheet1"))


x0 = Float64.(res_df.CD)


dev_stand0 = Statistics.std(x0)

merged_db = innerjoin(sampling_df, res_df, on = :SimCode)

merged_dbs = sort(merged_db,:CD)

Nsampl = size(sampling_df)[1]

plotly()
plot(xlabel = "Sample", ylabel="Parameter" )
plot!(1:1:Nsampl,merged_dbs.CD, markershape=:square, markersize=2, label="separation location [x/c]")
plot!(1:1:Nsampl,merged_dbs.Reat,markershape=:square, markersize=2, label="s1/10")

plot!(1:1:Nsampl,merged_dbs.s1 ./ 10,markershape=:square, markersize=2, label="s1/10")
# savefig("DU89_separation_s10_corr.pdf")


M = zeros(4,9)
corr_res = zeros(9)
for (i,corr_cmp) in enumerate([merged_dbs.CL,  merged_dbs.CD, merged_dbs.Sep,merged_dbs.Reat])


corr_res[1] = StatsBase.corspearman(Float64.(corr_cmp), Float64.(merged_dbs.U_freestream))
corr_res[2] =StatsBase.corspearman(Float64.(corr_cmp), Float64.(merged_dbs.TI))
corr_res[3] =StatsBase.corspearman(Float64.(corr_cmp), Float64.(merged_dbs.AoA))

corr_res[4] =StatsBase.corspearman(Float64.(corr_cmp), Float64.(merged_dbs.Viscosity_Ratio))
corr_res[5] =StatsBase.corspearman(Float64.(corr_cmp), Float64.(merged_dbs.σw1))
corr_res[6] =StatsBase.corspearman(Float64.(corr_cmp), Float64.(merged_dbs.α1))
corr_res[7] =StatsBase.corspearman(Float64.(corr_cmp), Float64.(merged_dbs.βstar))
corr_res[8] =StatsBase.corspearman(Float64.(corr_cmp), Float64.(merged_dbs.s1))
corr_res[9] =StatsBase.corspearman(Float64.(corr_cmp), Float64.(merged_dbs.C1))
M[i,:] = corr_res
end
var_names = ["U∞", "TI", "AoA", "μr", "σw1","α1","β*","s1", "C1"]
outputs = ["CL", "CD", "Separation", "Reattachment"]
heatmap(var_names,outputs, M, aspect_ratio=:equal, colormap=:turbo)

savefig("heatmap_uq_.pdf")



# plot!(1:1:110,merged_dbs.CD*100)

plot(ylabel="Cf", xlabel="x/c")
plot_cf(merged_dbs, 110;chord = 0.2, style=:dash, color=:blue, label="Reattach max")
plot_cf(merged_dbs, 1;chord = 0.2, style=:dash, color=:orange, label="Reattach min")
plot!(xlims=([0,1]),ylims=([-0.02,0.02]))

merged_dbs
gr()
plot(ylabel="Cf", xlabel="x/c")
plot_cf(merged_dbs, 1; chord = 0.2, style=:dash, color=:blue, label = "Reattachment point = 0.72c" )
plot_cf(merged_dbs, 110; chord = 0.2, style=:dash, color=:violet,label="Reattachment point = 0.95c" )

plot!(xlims=([0,1]),ylims=([-0.02,0.02]), legend=true)
savefig("DU89_UQ_cf_limits.pdf")


xd = Float64.(sampling_df.α1)

d = fit(Normal, xd)
xn = collect(LinRange(minimum(xd),maximum(xd),2000))
pd = pdf.(d, xn)
plot(xn,pd, label=false, linewidth=1.5)






# J_sample = Float64.(XLSX.readdata("Airfoil/SamplingN_$(airfoilname).xlsx", "Sheet1", "A2:I111")) #du89
J_sample = Float64.(XLSX.readdata("Airfoil/SamplingN_$(airfoilname).xlsx", "Sheet1", "A2:G73")) #sd7003

Ap1 = evaluate(J_sample,mop)'


yp = Ap1\x0


var_names

exp_val, dev_stand = pce_analysis(mop,yp)

dev_stand/dev_stand0


S = zeros(length(var_names))





#Sobol Indexes SD7003

sobol_df = Float64.((XLSX.readdata("Airfoil/Results/$(airfoilname)_results_Sobol.xlsx", "Sheet1","A2:C8")))

sobol_df

var_names = ["TI", "μr", "σw1","α1","β*", "s1", "C1"]

mn = [sobol_df[:,1]...,sobol_df[:,2]...,sobol_df[:,3]...] #[20, 35, 30, 35, 27,25, 32, 34, 20, 25]


sx = CategoricalArray(repeat(["E", "Separation", "Reattachment"], inner = 7))
levels!(sx, ["E", "Separation", "Reattachment"])

nam = CategoricalArray([var_names...,var_names...,var_names...])
levels!(nam, var_names)


splot = StatsPlots.groupedbar(nam, mn, group = sx, xlabel = "Variables",  ylabel = "Sobol Index",lw = 0)



savefig(splot,"Sobol_$(airfoilname).pdf")

xS = collect(0:1:length(var_names)-1)
bar_sob = bar(xS, S,xticks=(xS,var_names), label = false, ylabel="Sₜ")
savefig(bar_sob, "Sobol_Indexes_$(airfoilname).pdf")


plt_E = plot(ylabel="PDF", xlabel="E")
plot!(xn,pd, label=false, linewidth=1.5)

plot!(exp_val .* ones(10), LinRange(0,maximum(pd),10),label="expected value", linewidth=1.5)

plot!(exp_val .* ones(10) .+ 2*dev_stand, LinRange(0,maximum(pd),10),label="95% confidence interval", linewidth=1.1, linestyle=:dash, linecolor=:black)
plot!(exp_val .* ones(10) .- 2*dev_stand, LinRange(0,maximum(pd),10),label=false, linewidth=1.1, linestyle=:dash, linecolor=:black)




#Sobol Indexes DU89

sobol_df = Float64.((XLSX.readdata("Airfoil/Results/$(airfoilname)_results_Sobol.xlsx", "Sheet1","C2:E10")))

var_names = ["U_freestream", "TI", "AoA", "Viscosity_Ratio", "σw1","α1","βstar","s1", "C1"]

mn = [sobol_df[:,1]...,sobol_df[:,2]...,sobol_df[:,3]...] #[20, 35, 30, 35, 27,25, 32, 34, 20, 25]


sx = CategoricalArray(repeat(["E", "Separation", "Reattachment"], inner = length(var_names)))
levels!(sx, ["E", "Separation", "Reattachment"])

nam = CategoricalArray([var_names...,var_names...,var_names...,var_names...,var_names...])
levels!(nam, var_names)


splot = StatsPlots.groupedbar(nam, mn, group = sx, xlabel = "Variables",  ylabel = "Sobol Index",lw = 0)



savefig(splot,"Sobol_$(airfoilname).pdf")


d = fit(Normal, x)
xn = collect(LinRange(minimum(x),maximum(x),2000))
pd = pdf.(d, xn)

plt_E = plot(ylabel="PDF", xlabel="E")
plot!(xn,pd, label=false, linewidth=1.5)

plot!(exp_val .* ones(10), LinRange(0,maximum(pd),10),label="expected value", linewidth=1.5)

plot!(exp_val .* ones(10) .+ 2*dev_stand, LinRange(0,maximum(pd),10),label="95% confidence interval", linewidth=1.1, linestyle=:dash, linecolor=:black)
plot!(exp_val .* ones(10) .- 2*dev_stand, LinRange(0,maximum(pd),10),label=false, linewidth=1.1, linestyle=:dash, linecolor=:black)