using Distributions, XLSX, DataFrames, Distributions, StatsPlots, CategoricalArrays, JLD2
using Plots.Measures

include(joinpath(@__DIR__,"..","..","UncertaintyUtils.jl"))
include(joinpath(@__DIR__, "UQAirfoilPost.jl"))


dirname = joinpath(@__DIR__,"..","Results_2024","Results_ParisOspital","Results")
dirnameCFResults = joinpath(@__DIR__,"..","Results_2024","Results_ParisOspital","ResultsCF")

file_sampling = joinpath(@__DIR__,"..","Results_2024","Sampling_DU89.XLSX")
file_samplingN =  joinpath(@__DIR__,"..","Results_2024","SamplingN_DU89.XLSX")

mop = load("MOP_DU89.jld2")["mop"]



uqairfoil = UQAirfoil(dirname, dirnameCFResults; chord = 0.2, confidence=0.05, writeCF=false, compute_outliers=true)
uq_cumulative = CumulativeUQAirfoil(uqairfoil)



gr()

flabel= ("CL", "CD", "separation [x/c]", "reattachment [x/c]")
fieldnames(typeof(uq_cumulative))

gr()
default(linewidth=2)
PLT_CUMULATIVE = Plots.Plot[]
for (i,fname) in enumerate(fieldnames(typeof(uq_cumulative)))
    field_cumulate = getfield(uq_cumulative,fname)
    label_expcted = false
    label_confidence = false
    if i == 1
        label_expcted = "expected value"
        label_confidence = "95% confidence interval"
    end

    pltf = plot(field_cumulate.mean, linecolor=:deepskyblue3, linestyle=:solid, label = label_expcted)
    plot!(field_cumulate.confidence_interval_top, linecolor=:black, linestyle=:dash, label = label_confidence)
    plot!(field_cumulate.confidence_interval_bottom, linecolor=:black, linestyle=:dash, label= false)
    plot!(xlabel= "number of simulations", ylabel=flabel[i])
    push!(PLT_CUMULATIVE,pltf)
end

plot(PLT_CUMULATIVE..., layout=(2,2), size=(1000,700), left_margin = 3mm)


savefig("Cumulative_convergence.pdf")


plt_CL = plot_distibution(uqairfoil.CL.all_values, uqairfoil.CL.confidence_interval, "CL"; labelling=true)
plt_CD = plot_distibution(uqairfoil.CD.all_values, uqairfoil.CD.confidence_interval, "CD")

plt_sep = plot_distibution(uqairfoil.xsep.all_values, uqairfoil.xsep.confidence_interval, "separation [x/c]")
vline!([0.59], label="experimental value", linewidth=1.6, linestyle=:solid, linecolor=:black)
vspan!([0.57,0.62]; color=:grey, alpha=0.3, label="experimental error band")

plt_reatt = plot_distibution(uqairfoil.xreatt.all_values, uqairfoil.xreatt.confidence_interval, "reattachment [x/c]")
vline!([0.71], label=false, linewidth=1.6, linestyle=:solid, linecolor=:black)
vspan!([0.68,0.75]; color=:grey, alpha=0.3, label=false)

plot(plt_CL,plt_CD,plt_sep,plt_reatt, layout=(2,2), size=(1000,700), left_margin = 3mm)

savefig("Distributions.pdf")



##### CF uncertainty

cfx, cf_mean,cf_lower_perc, cf_top_per = compute_CF_UQ(uqairfoil::UQAirfoil;   confidence=0.05, Nx = 250 )

gr()
plot(cfx,cf_mean, linecolor=:deepskyblue3, linestyle=:solid, label = "expected value" )
plot!(cfx, cf_lower_perc, fillrange = cf_top_per, fillalpha = 0.55, linewidth=0.0, c = 1, label = "95% confidence interval", legend = :topleft)
plot!(ylims=([-0.01,0.02]), xlabel="x/c")
Plots.savefig("CF_UQ.pdf")

###PARALLEL PLOTS
uqvar_fname = joinpath(@__DIR__,"..","Results_2024","Sampling_DU89.xlsx")
uq_vars = UQVariables(uqvar_fname)
uqairfoil_all = UQAirfoil(dirname, dirnameCFResults; chord = 0.2, confidence=0.05, writeCF=false, compute_outliers=true)

df_IDX = zeros(Int64,length(uqairfoil_all.SimCode[uqairfoil.not_outliers_idx]))

for (i,scode) in enumerate(uqairfoil_all.SimCode[uqairfoil.not_outliers_idx])
    println(i)
    df_IDX[i] = findall(x->x== scode, uq_vars.df.SimCode)[1]
end


using PlotlyJS, DataFrames, PlotlyBase
using Plotly

plotly()

df = hcat(uq_vars.df[[df_IDX...],:],DataFrame(:CL=>uqairfoil_all.CL.all_values),DataFrame(:CD=>uqairfoil_all.CD.all_values),
DataFrame(:xsep=>uqairfoil_all.xsep.all_values),DataFrame(:xreatt=>uqairfoil_all.xreatt.all_values)) 



in_out_names = ["U₀", "TI", "AoA", "μr", "σw1","α1","βstar","s1", "C1", "CL", "CD", "sep", "reatt"]


trace_inout= parcoords(;line = attr(color=df.CD, colorscale=:turbo),
dimensions = [attr(label=in_out_names[i], values=val,     

tickvals = [minimum(val), maximum(val)],    ticktext = [string(round(minimum(val);sigdigits =3)), string(round(maximum(val);sigdigits =3))])

for (i,val) in enumerate(eachcol(df[:,[1:9...,11:end...]]))] )
plot_parallel_in_out = PlotlyJS.plot(trace_inout, Layout(font = attr(size = 20) ))

PlotlyJS.savefig(plot_parallel_in_out, "plot_parallel_in_out.html")


trace_out = parcoords(;line = attr(color=df.CD, colorscale=:turbo),
dimensions = [attr(label = "CL", values = df.CL),attr(label = "CD", values = df.CD), 
attr(label = "xsep", values = df.xsep),attr(label = "xreatt", values = df.xreatt),])               
plot_parallel_out = PlotlyJS.plot(trace_out, )

PlotlyJS.savefig(plot_parallel_out, "plot_parallel_out.html")


### POLYNOMIAL CHAOS EXPANSION
J_sample = Matrix(UQVariables(file_sampling)(uqairfoil.SimCode[uqairfoil.not_outliers_idx])[:,1:end-1])
J_sample_eval = Matrix(UQVariables(file_samplingN)(uqairfoil.SimCode[uqairfoil.not_outliers_idx])[:,1:end-1])

Ap1 = evaluate(J_sample_eval,mop)'


yp = Ap1\uqairfoil.xreatt.all_values
exp_value, dev_std = pce_analysis(mop,yp)



####SOBOL

gr()

#Sobol Indexes
airfoilname ="DU89"
sobol_df = Float64.((XLSX.readdata( joinpath(@__DIR__,"..","Results_2024", "DU89_results_Sobol.xlsx"), "Sheet1","A2:D10")))
var_names = ["U₀", "TI", "AoA", "μr", "σw1","α1","βstar","s1", "C1"]
num_vars = length(var_names) #number of variables





## Create Histogram Sobol Indexes
mn = [sobol_df[:,1]...,sobol_df[:,2]...,sobol_df[:,3]...,sobol_df[:,4]...] #[20, 35, 30, 35, 27,25, 32, 34, 20, 25]

sx = CategoricalArray(repeat(["CL","CD", "Separation", "Reattachment"], inner = num_vars))
levels!(sx, ["CL","CD", "Separation", "Reattachment"])

nam = CategoricalArray([var_names...,var_names...,var_names...,var_names...])
levels!(nam, var_names)


splot = StatsPlots.groupedbar(nam, mn, group = sx, xlabel = "Variables",  ylabel = "Sobol Index",lw = 0)

savefig(splot,"Sobol_DU89.pdf")




#### Comparison MC vs PCE

## Cumulative PCE
include(joinpath(@__DIR__, "PCE_CUMULATIVE.jl"))

CL_pce_cumulative = PCECumulative(Ap1,uqairfoil.CL.all_values,mop)
CD_pce_cumulative = PCECumulative(Ap1,uqairfoil.CD.all_values,mop)
Xsep_pce_cumulative = PCECumulative(Ap1,uqairfoil.xsep.all_values,mop)
Xreatt_pce_cumulative = PCECumulative(Ap1,uqairfoil.xreatt.all_values,mop)

PCE_cumulative_results = [CL_pce_cumulative, CD_pce_cumulative, Xsep_pce_cumulative, Xreatt_pce_cumulative]

gr()

PLT_CUMULATIVE_PCE_MC = Plots.Plot[]
for (i,fname) in enumerate(fieldnames(typeof(uq_cumulative)))
    field_cumulate = getfield(uq_cumulative,fname)
    label_expcted = false
    label_confidence = false
    label_expcted_pce = false
    label_confidence_pce = false

    if i == 1
        label_expcted = "expected value MC"
        label_confidence = "95% confidence interval MC"
    elseif i ==2
           label_expcted_pce = "expected value PCE"
        label_confidence_pce = "95% confidence interval PCE"
    end

    pce_res = PCE_cumulative_results[i]

    pltf = plot(field_cumulate.mean, linecolor=:blue, linestyle=:solid, label = label_expcted)
    plot!(field_cumulate.confidence_interval_top, linecolor=:black, linestyle=:dash, label = label_confidence)
    plot!(field_cumulate.confidence_interval_bottom, linecolor=:black, linestyle=:dash, label= false)
    

    plot!(pce_res.mean, linecolor=:chartreuse4, linestyle=:solid, label = label_expcted_pce)
    plot!(pce_res.confidence_interval_top, linecolor=:green1, linestyle=:dash, label = label_confidence_pce)
    plot!(pce_res.confidence_interval_bottom, linecolor=:green1, linestyle=:dash, label= false)
    
    
    
    plot!(xlabel= "number of simulations", ylabel=flabel[i], ylims=([0.9*minimum(field_cumulate.confidence_interval_bottom),1.1*maximum(field_cumulate.confidence_interval_top)]))



    push!(PLT_CUMULATIVE_PCE_MC,pltf)
end

plot(PLT_CUMULATIVE_PCE_MC..., layout=(2,2), size=(1000,700), left_margin = 3mm)

savefig("PCE_vs_MC_DU89.pdf")

