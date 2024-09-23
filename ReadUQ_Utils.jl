using Statistics, CSV, DataFrames, XLSX, KissSmoothing





function plot_cf(sampling_df::DataFrame, idx::Int64;chord = 1.0, style=:solid, color=:blue, label=false)
    cname = sampling_df.SimCode[idx]
    fname = joinpath("Results1","CFSim" * cname * ".csv")
    curve = read_CF_file(fname)
    cf_clean = denoise(curve.cf, factor = 0.5)
    plot!(curve.x ./ chord,cf_clean; linestyle=style, linecolor=color, label = label, linewidth=1.5)
end


function read_CF_file(filename)
data = CSV.read(filename, DataFrame, skipto = 2, header=["x","cf"])
curve = (DataFrame(Float64.(data[2:end,:])))
curve = filter(d->d.cf>-50,curve)
curve = sort(curve, "x")
return curve
end


function find_lsb_points(filename; chord = 1.0)

curve = read_CF_file(filename)

s = sign.(curve.cf[2:end] .* curve.cf[1:end-1])  .== -1

idxs = findall(v-> v == true,s)

xv = curve.x[idxs] ./ chord

id_sep = findfirst(x->x>0.3,xv)
id_reat = findlast(x->x>0.3,xv)

return xv[id_sep], xv[id_reat]
end

function get_CL(filename)
    data = CSV.read(filename, DataFrame, skipto = 2, header=["iter","CL"])
    CL = Statistics.mean(data.CL[end-100,end])
    return CL
end

function get_CD(filename)
    data = CSV.read(filename, DataFrame, skipto = 2, header=["iter","CD"])
    CD = Statistics.mean(data.CD[end-100,end])
    return CD
end



function read_UQ_results(filename::String; chord = 1.0)

SamplingFile = DataFrame(XLSX.readtable(filename, "Sheet1"))

Sep = Float64[]
Reat = Float64[]
CL = Float64[]
CD = Float64[]

for cname in String.(SamplingFile.SimCode)
    sep_tmp, reat_tmp = find_lsb_points(joinpath("Results1","CFSim$cname.csv"); chord = chord)
    CL_tmp = get_CL(joinpath("Results1","CLSim$cname.csv"))
    CD_tmp = get_CD(joinpath("Results1","CDSim$cname.csv"))
    push!(Sep,sep_tmp)
    push!(Reat,reat_tmp)
    push!(CL,CL_tmp)
    push!(CD,CD_tmp)

end

res_df = DataFrame(:Sep=>Sep, :Reat=> Reat,:CL=>CL, :CD=>CD)
return res_df
end

