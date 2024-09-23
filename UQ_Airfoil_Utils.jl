"""
return J_sample_eval, J_sample. J_sample is the random sampling in having μ=0.0, σ=1.0, , while J_sample_eval is the corresponding sampling in the physical specified domain.

Internally it verify that the values respect the following condition
σ*3-μ < val  <σ*3+μ

In the idx there are the gaussian distributions

And apply the half-normal distribution at the variable at index half_normal
"""
function get_Samples(mop, Ns::Int64, multi_pce,multi_dist, multi_mean, multi_σ; idxs=[1,2,3,5,6,7], half_normal=6)
    nx = length(multi_mean)
    J_sample = zeros(Ns,nx)
    J_sample_eval = zeros(Ns,nx)
    
    nsample = 0
    J_sample_tmp = PolyChaos.sampleMeasure(Ns*5,mop) #*4 to be sure enough data after discarding not accepted values
    J_sample_tmp_eval = zeros(size(J_sample_tmp))


 

idx_cond = Vector{Float64}[]

for (i,(mp,md)) in enumerate(zip(multi_pce,multi_dist))
    J_sample_tmp_eval[:,i] = evaluatePCE(mp, J_sample_tmp[:,i],md)
end



for i in idxs
        idx_tmp = abs.(evaluatePCE(multi_pce[i], J_sample_tmp[:,i], multi_dist[i]) .- multi_mean[i] ) .< 3 * multi_σ[i]
    if i == half_normal
        idx_tmp2 = evaluatePCE(multi_pce[i], J_sample_tmp[:,i], multi_dist[i]) .>= multi_mean[i] 
        idx_tmp = idx_tmp .* idx_tmp2
    end
    push!(idx_cond, idx_tmp)
end

cond_true = Float64[]
for i = 1:1:length(idx_cond[1])
    push!(cond_true,prod(map(e-> e[i], idx_cond)))
end


nsample = 0
inner_iter = 1
while nsample < Ns
    

    if cond_true[inner_iter] >0
        println(nsample)
        nsample = nsample+1
        J_sample[nsample,:] = J_sample_tmp[inner_iter,:]
        J_sample_eval[nsample,:] = J_sample_tmp_eval[inner_iter,:]

    end
    inner_iter = inner_iter +1
end

return J_sample_eval, J_sample

end


function export_samples(fname::String, J_sample_eval::Matrix, var_names::Vector)
    J_sample_eval_df = DataFrame( var_names[1] => J_sample_eval[:,1],var_names[2] => J_sample_eval[:,2],var_names[3] => J_sample_eval[:,3],
    var_names[4] => J_sample_eval[:,4],var_names[5] => J_sample_eval[:,5],var_names[6] => J_sample_eval[:,6], var_names[7] => J_sample_eval[:,7], var_names[8] => J_sample_eval[:,8],var_names[9] => J_sample_eval[:,9])
    XLSX.writetable(fname, J_sample_eval_df)
    return fname
end


