using Random
using Statistics


alphaa = 0.05        # Significance level for 95% confidence interval

# Generate samples


samples = CDres
# Estimate the mean
mean_estimate = mean(samples)

# Compute the confidence interval

lower_percentile = 100 * (alphaa / 2)
upper_percentile = 100 * (1 - alphaa / 2)

confidence_interval = quantile(samples, [lower_percentile / 100, upper_percentile / 100])

confidence_interval
println("Estimated Mean: $mean_estimate")
println("95% Confidence Interval: $confidence_interval")


nbin=10

histogram(samples) 

vline!(mean_estimate .* ones(100), label="expected value", linewidth=1.5)

vline!(confidence_interval[1] .* ones(100), label="95% confidence interval", linewidth=1.1, linestyle=:dash, linecolor=:black)
vline!(confidence_interval[2] .* ones(100), label=false, linewidth=1.1, linestyle=:dash, linecolor=:black)

