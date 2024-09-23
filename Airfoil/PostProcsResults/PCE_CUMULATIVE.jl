using Statistics
include(joinpath(@__DIR__, "UQAirfoilPost.jl"))

function PCECumulative(AP,out_val,mop)
    Exp_Value = zeros(length(out_val))
    Dev_Std = zeros(length(out_val))
    confidence_interval_bottom = zeros(length(out_val))
    confidence_interval_top = zeros(length(out_val))

    for i in eachindex(out_val)
        Ap1 = AP[1:i,:]
        yp = Ap1\out_val[1:i]
        Exp_Value[i], Dev_Std[i] = pce_analysis(mop,yp)
    end
    confidence_interval_bottom = Exp_Value - 2*Dev_Std
    confidence_interval_top = Exp_Value + 2*Dev_Std
    CumulativeStatisticalOutput(Exp_Value, confidence_interval_bottom, confidence_interval_top)
end