using TuringGLM, BenchmarkTools

using TuringGLM, BenchmarkTools, DataFrames, Turing, Plots
# Generate some sample data
n = 1000
x1 = rand(n)
x2 = rand(n)
x3 = rand(n)
y = 3 .+ 2 .* x1 .- 4 .* x2 + 0.5 .* randn(n) .+ 0.5 .* x3
df = DataFrame(y=y, x1=x1, x2=x2, x3=x3)
# Define the formula for linear regression
fm_linear = @formula(y ~ x1 + x2 + x3)

# Instantiate the Turing model for linear regression
model_linear = turing_model(fm_linear, df)

# Number of samples for NUTS
n_samples = 2000

# Benchmarking the sampling for linear model
linear_bench = @benchmark begin
    chns_linear = sample(model_linear, NUTS(), $n_samples)
end
# Define the formula for robust regression (assuming same formula)
fm_robust = @formula(y ~ x1 + x2 + x3)

# Instantiate the Turing model for robust regression (may require specific configurations)
model_robust = turing_model(fm_robust, df, model=TDist)

# Benchmarking the sampling for robust model
robust_bench = @benchmark begin
    chns_robust = sample(model_robust, NUTS(), $n_samples)
end
# Extracting and displaying the time taken
linear_time = median(linear_bench.times) / 1e9  # Converting to seconds
robust_time = median(robust_bench.times) / 1e9  # Converting to seconds

println("Time taken for Linear Model: ", linear_time, " seconds")
println("Time taken for Robust Model: ", robust_time, " seconds")
println("Ratio of times: ", robust_time / linear_time)
