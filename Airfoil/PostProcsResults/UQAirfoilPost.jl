
    

#### POST PROCESS UQ CSV RESULTS ####### 


using Statistics, CSV, DataFrames, Statistics, XLSX
using Plots
using Distributions
using Interpolations

function convert_num2str(n::Int64)
    if n < 1 || n >= 26^3
      error("Input must be between 1 and 26^3 - 1 (inclusive).")
    end
  
    # Initialize an array for the letter sequence
    letters = ['A', 'A', 'A']
  
    # Convert to 0-based index
    n -= 1
  
    for i in 1:3
        letters[4 - i] = Char('A' + (n % 26))
        n ÷= 26
    end
  
    return join(letters)
end

function write_CF_file_curve(filename, dirname, dirnameCFResults; chord = 0.2)
    filenameread = joinpath(dirname, filename)
    data = CSV.read(filenameread, DataFrame, skipto = 2, header=["x","cf"])[:,1:2]
    curve = filter(d->d.cf>-50,data)
    curve = sort(curve, "x")

    curve.x = curve.x ./ chord
    outputfile = joinpath(dirnameCFResults, filename)
    println("Writing $filename")
    CSV.write(outputfile, curve)
end



function find_separation_reattachment_point(filename::String)
    data = CSV.read(filename, DataFrame)


    s = sign.(data.cf[2:end] .* data.cf[1:end-1])  .== -1

    idxs = findall(v-> v == true,s)

    data.x[idxs].>0.25

    id = findall(x->x>0.25, data.x[idxs])
    sep = 0.0
    reatt = 1.0

    if length(id)>= 2 
        sep = data.x[idxs[id[1]]]
        reatt = data.x[idxs[id[2]]]
    elseif length(id) == 1
        sep = data.x[idxs[id[1]]]
    end

    return data, sep, reatt
end

function find_separation_reattachment_point(filenames::Vector{String})
    n = length(filenames)
    xsep = zeros(n)
    xreatt = zeros(n)

    for (i,filename) in enumerate(filenames)
        xsep[i],xreatt[i]= find_separation_reattachment_point(filename)
    end
    return xsep,xreatt
end

function find_outilers(xreatt::Vector{Float64}, compute_outliers; threshold=1.0 )
    idx_out = []
    idx_n_out = collect(eachindex(xreatt))

    if compute_outliers
        idx_out = findall(x->x >= threshold, xreatt)
        idx_n_out = findall(x->x < threshold, xreatt)
    end

    return idx_n_out, idx_out
end

function get_CL(filename)
    data = CSV.read(filename, DataFrame, skipto = 2, header=["iter","CL"])[:,1:2]
    return Statistics.mean(data.CL[end-100:end])
end


function get_CD(filename)
    data = CSV.read(filename, DataFrame, skipto = 2, header=["iter","CD"])[:,1:2]
    return Statistics.mean(data.CD[end-100:end])
end



function plot_distibution(samples, confidence_interval, xlabel; labelling = false)
    distributions_univarate = [Cauchy{Float32}, Normal, Uniform, Gamma, LogNormal]
    fd = Distribution[]
    for dist in distributions_univarate
        try
            fdtmp = fit.(dist, Ref(samples))
            push!(fd, fdtmp)
        catch
        end
    end

       label_expected = false
        label_confidence =  false
    if labelling
        label_expected = "expected value"
        label_confidence =  "95% confidence interval"
    end

    d = fd[findmax(loglikelihood.(fd, Ref(samples)))[2]]
    println("Distribution = $(typeof(d)), $d")
    xpdf = collect( LinRange(minimum(samples).* 0.95 ,maximum(samples).*1.05,1000))

    histogram(samples, normalize=:pdf, label=false, linecolor=:deepskyblue3,color=:deepskyblue3) 
    plot!(xpdf, pdf.(d, xpdf), linewidth=2,linecolor=:red,label=false)
    vline!(mean(samples) .* ones(100), label=label_expected, linewidth=1.1, linestyle=:solid, linecolor=:blue)
    vline!(confidence_interval[1] .* ones(100), label=label_confidence, linewidth=1.1, linestyle=:dash, linecolor=:blue)
    vline!(confidence_interval[2] .* ones(100), label=false, linewidth=1.1, linestyle=:dash, linecolor=:blue)

    plot!(xlabel=xlabel, ylabel="P(X)")

end



struct StatisticalOutput
    expected_value::Float64 ## mean
    confidence_interval::Vector ## 95 confidence interval
    all_values::Vector    

end



function StatisticalOutput(samples::Vector, confidence, idx_n_out)
    samples=samples[idx_n_out]
    lower_percentile = 100 * (confidence / 2)
    upper_percentile = 100 * (1 - confidence / 2)
    confidence_interval = quantile(samples, [lower_percentile / 100, upper_percentile / 100])
    StatisticalOutput(Statistics.mean(samples),confidence_interval, samples)
end


struct UQAirfoil
    CL::StatisticalOutput
    CD::StatisticalOutput
    xsep::StatisticalOutput
    xreatt::StatisticalOutput
    SimCode::Vector
    CF::Vector
    not_outliers_idx::Vector
    outliers_idx::Vector
end

struct SubUQAirfoil
    CL::Float64
    CD::Float64
    xsep::Float64
    xreatt::Float64
    SimCode::String
    CF::DataFrame
end



function find_simcode(dirname)
    idx_start =findfirst.("Sim",readdir(dirname))[1][end]
    unique_strings = unique(map(x->x[idx_start+1:idx_start+3], readdir(dirname)))
end

function UQAirfoil(dirname, dirnameCFResults; chord = 0.2, confidence=0.05, writeCF=true, compute_outliers=false)
    files_results_names = readdir(dirname)
    files_CF0 = files_results_names[findall(occursin.("CFSim",files_results_names))]
    writeCF && write_CF_file_curve.(files_CF0, dirname, dirnameCFResults; chord = chord)
    files_CF = readdir(dirnameCFResults)

    unique_strings = find_simcode(dirname)
    files_CF = readdir(dirnameCFResults)
    
    CL,CD,xsep,xreatt= [zeros(length(unique_strings)) for i = 1:4]
    CF = DataFrame[]
    
    for (i,scode) in enumerate(find_simcode(dirname))
        fileCL = joinpath(dirname, files_results_names[findall(occursin.("CLSim$scode",files_results_names))][1])
        fileCD = joinpath(dirname,files_results_names[findall(occursin.("CDSim$scode",files_results_names))][1])
        fileCF = joinpath(dirnameCFResults ,files_CF[findall(occursin.("CFSim$scode",files_CF))][1])
        data, xsep[i], xreatt[i] = find_separation_reattachment_point(fileCF)
        push!(CF,data)
        CL[i] = get_CL(fileCL)
        CD[i] = get_CD(fileCD)
    end

    idx_n_out, idx_out = find_outilers(xreatt,compute_outliers; threshold=1.0)
    UQAirfoil(
    StatisticalOutput(CL, confidence, idx_n_out),
    StatisticalOutput(CD, confidence, idx_n_out),
    StatisticalOutput(xsep, confidence, idx_n_out),
    StatisticalOutput(xreatt, confidence, idx_n_out),
    unique_strings,
    CF,idx_n_out,idx_out)
    
end

function (uqairfoil::UQAirfoil)(xs::Vector{String})
    idx_xs = [findall(x -> x == s, uqairfoil.SimCode)[1] for s in xs]
   
    subuq = SubUQAirfoil[]
    for ii in idx_xs

        subuqi = SubUQAirfoil(
        uqairfoil.CL.all_values[ii],
        uqairfoil.CD.all_values[ii],
        uqairfoil.xsep.all_values[ii],
        uqairfoil.xreatt.all_values[ii],
        uqairfoil.SimCode[ii],
        uqairfoil.CF[ii])

        push!(subuq, subuqi)
    end

    return subuq
end

### CF uncertainty
function compute_CF_UQ(uqairfoil::UQAirfoil;   confidence=0.05, Nx = 250 )
    xx = collect(LinRange(0.0,1.0,Nx))
    CF = zeros(Nx, length(uqairfoil.CF))

    for (i,df_cf) in enumerate(uqairfoil.CF)
        li = linear_interpolation(df_cf.x,df_cf.cf,extrapolation_bc=Line())
        CF[:,i] = li.(xx)
    end

    CFxx = [StatisticalOutput(CF[i,:], confidence, uqairfoil.not_outliers_idx) for i = 1:Nx]
    return xx,map(so-> so.expected_value, CFxx),map(so-> so.confidence_interval[1], CFxx),map(so-> so.confidence_interval[2], CFxx)
end

struct CumulativeStatisticalOutput
    mean::Vector
    confidence_interval_bottom::Vector
    confidence_interval_top::Vector

end


function CumulativeStatisticalOutput(so::Vector{StatisticalOutput})
    CumulativeStatisticalOutput(map(x->x.expected_value, so),
    map(x->x.confidence_interval[1], so),
    map(x->x.confidence_interval[2], so))
end

struct CumulativeUQAirfoil
    CL::CumulativeStatisticalOutput
    CD::CumulativeStatisticalOutput
    xsep::CumulativeStatisticalOutput
    xreatt::CumulativeStatisticalOutput

end






function CumulativeUQAirfoil(uqairfoil::UQAirfoil)
    CL_SO = StatisticalOutput[]
    CD_SO = StatisticalOutput[]
    xsep_SO = StatisticalOutput[]
    xreatt_SO = StatisticalOutput[]
    
    confidence=0.05
    idx_n_out = uqairfoil.not_outliers_idx
    for IDX in eachindex(uqairfoil.CL.all_values)
        idxs = collect(1:IDX)
        push!(CL_SO, StatisticalOutput(uqairfoil.CL.all_values[idxs], confidence, intersect(idxs,idx_n_out)))
        push!(CD_SO,StatisticalOutput(uqairfoil.CD.all_values[idxs], confidence, intersect(idxs,idx_n_out)))
        push!(xsep_SO,StatisticalOutput(uqairfoil.xsep.all_values[idxs], confidence, intersect(idxs,idx_n_out)))
        push!(xreatt_SO,StatisticalOutput(uqairfoil.xreatt.all_values[idxs], confidence, intersect(idxs,idx_n_out)))
    end
    CumulativeUQAirfoil(
    CumulativeStatisticalOutput(CL_SO), CumulativeStatisticalOutput(CD_SO), 
    CumulativeStatisticalOutput(xsep_SO), CumulativeStatisticalOutput(xreatt_SO))
end


struct UQVariables
    df::DataFrame
end


function UQVariables(filename)
    xf = DataFrame(XLSX.readtable(filename, "Sheet1"))
    nsim = size(xf)[1]
    xfnew = Float64.(xf[:,1:end-1])
    xfnew[!, :SimCode] =  String.(xf[:,end])
    return UQVariables(xfnew)
end

function (uqvar::UQVariables)(s::String)
    return filter(r->r.SimCode==s, uqvar.df)
end



function (uqvar::UQVariables)(s::Vector{String})
    df = uqvar(s[1])
    for si in s[2:end]
        df_tmp =  uqvar(si)
        df = vcat(df,df_tmp)
    end
    return df
end


# end





