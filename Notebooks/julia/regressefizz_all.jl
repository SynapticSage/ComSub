using DataFrames, Statistics, LinearAlgebra, Serialization
using Interpolations, Infiltrator
using MAT
using Plots
unicodeplots()
using Distributed
using ProgressMeter
using DataStructures: OrderedDict
using StatsPlots
using PyCall
@pyimport sklearn.metrics as skm
using SearchSortedNearest
include("/Volumes/MATLAB-Drive/Shared/Notebooks/julia/regess_funcs.jl")
folder = "/Volumes/MATLAB-Drive/Shared/figures/midpattern=true/data/"
file = "/Volumes/MATLAB-Drive/Shared/figures/midpattern=true/data/$(animal)_coherence_true.mat";
overwrite=false;
trades = [r"phi_cos" => "phi-cos", r"phi_sin" => "phi-sin",
    r"low_gamma" => "low-gamma", r"mid_gamma" => "mid-gamma", 
    r"high_gamma" => "high-gamma", r"wpli_avg" => "wpli-avg"]


usegpu = false;
# n_chains = 16;
n_chains = 1;
# Number of samples
n_samples = Int(round(30/n_chains));
numsamprow = 25_000;
# Threads.nthreads() = n_chains;
# addprocs(n_chains-1)
@everywhere using TuringGLM, MCMCChains
if usegpu
    @everywhere using CUDA
    @everywhere using DynamicPPL
end
println("n_chains: $n_chains")
println("n_samples: $n_samples")
println("numsamprow: $numsamprow")

# Define frequency bands globally
bands = OrderedDict(
    "delta"      => (0.5, 4.0),
    "theta"      => (6.0, 12.0),
    "beta"       => (13.0, 30.0),
    "low_gamma"  => (30.0, 50.0),
    "mid_gamma"  => (50.0, 80.0),
    "high_gamma" => (80.0, 100.0),
    "epsilon"    => (100.0, 120.0),
    "ripple"     => (150.0, 250.0)
)

"""
    predict_and_compute_stats(model, test_data, predictors, all_predictors,
			      regressed)

Predicts on the test data using the model and computes summary statistics
such as accuracy, R2, MSE, precision, recall, and F1 score.

# Arguments
- `model`: The model to use for prediction
- `test_data`: The test data to predict on
- `predictors`: The predictors to use for prediction
- `all_predictors`: All the predictors in the dataset
- `regressed`: The column to regress on
- `paramrange`: The range of parameters to use for prediction
# Returns
- `accuracy`: The accuracy of the predictions
- `r2`: The R2 score of the predictions
- `mse`: The mean squared error of the predictions
- `predictions`: The predictions
- `prec`: The precision of the predictions
"""
function predict_and_compute_stats_v2(model,
    test_data, predictors, all_predictors, regressed)
    print("Predicting on test data...")
    paramrange = findall(in.(predictors, [all_predictors]))
    # Extract relevant coefficients
    α = model[:α]  # Intercept
    nms = [n for n in names(model) if contains(string(n), "β")]
    β = model[nms].value[:,:,1]  # Coefficients
    if paramrange !== nothing
	    β = β[:, paramrange]
    end
    # Predict on test data
    X_test = Matrix(test_data[:, predictors])
    # Ensure that the dimensions are consistent
    try
	@assert(!isempty(predictors))
@assert(size(X_test, 2) == size(β, 2))
	@assert(size(X_test, 1) == size(test_data, 1))
    catch
	@infiltrate
    end
    predictions = zeros(size(X_test,1), size(α,1))
    for i in eachindex(α)
	    predictions[:,i] = α[i] .+ (X_test * β[i, :].data)
    end
    # Compute summary stats
    true_values = test_data[:, regressed]  # Assuming "Behavior" is your target column
    if length(unique(true_values)) == 2
	    predictions = logistic.(predictions)
    end
    # R-squared
    residuals = true_values .- predictions
    residuals = residuals .^ 2;
    ss_res = vec(sum(residuals, dims=1))
    ts = (true_values .- mean(true_values)) .^ 2
    ts = repeat(ts, 1, size(predictions,2))
    ss_tot = vec(sum(ts, dims=1))
    r2 = 1 .- (ss_res ./ ss_tot)
    # Accuracy
    binary_preds = Int64.(round.(predictions))
    true_values = ones(Int64,size(binary_preds)) .* Int64.(round.(true_values));
    accuracy = vec(mean(binary_preds .== true_values, dims=1))
    bp = binary_preds[:,1]
    tv = true_values[:,1]
    prec = zeros(size(binary_preds,2))
    rec = zeros(size(binary_preds,2))
    @showprogress "prec and recall" for (i, (bp, tv)) in enumerate(zip(eachcol(binary_preds), eachcol(true_values)))
	    prec[i] = skm.precision_score(bp, tv, average="weighted")
	    rec[i] = skm.recall_score(bp, tv, average="weighted")
    end
    f1_score = 2 .* ((prec .* rec) ./ (prec.+ rec))
    # binarize = tv->skp.label_binarize(reshape(tv,length(tv),1), classes=unique(tv))
    # roc_auc = [skm.roc_auc_score(binarize(bp), binarize(tv), average="weighted") for (bp, tv) in zip(eachcol(binary_preds), eachcol(true_values))]
    # MSE
    mse = vec(mean(residuals .^ 2, dims=1))
    print("Accuracy: ", mean(accuracy))
    print("R2: ", mean(r2))
    print("MSE: ", mean(mse))
    return (;accuracy, r2, mse, predictions, prec, rec, f1_score)
end

function generate_frequency_bands(frequencies::Vector{Float64}, downsampling::Int=1)
    # Initialize an empty OrderedDict to hold the new bands
    new_bands = OrderedDict()
    
    # Iterate over each original band
    for (band_name, (lower, upper)) in bands
        # Filter the frequencies that fall within the current band
        filtered_frequencies = filter(f -> lower <= f <= upper, frequencies)
        
        # Apply downsampling
        for i in 1:downsampling:length(filtered_frequencies)
            start_idx = i
            end_idx = min(i + downsampling - 1, length(filtered_frequencies))
            
            new_band_name = string(band_name, "__", div(i - 1, downsampling) + 1)
            new_bands[new_band_name] = (filtered_frequencies[start_idx], filtered_frequencies[end_idx])
        end
    end
    
    return new_bands
end




function compute_band_means(efizza_data, frequencies)
    band_means = Dict()
    for (band_name, (low, high)) in bands
        idx = findall(x -> x >= low && x <= high, frequencies)
        band_means[band_name] = mean(efizza_data[:, idx], dims=2)
    end
    return band_means
end

function compute_band_means_for_efizz(efizz::Dict{String, Any})
    # Extract frequencies and time from efizz
    frequencies = efizz["f"][:]
    time_length = length(efizz["t"][:])

    # Initiate a new dictionary to store the band means
    band_means_dict = Dict()
    bands = generate_frequency_bands(frequencies, 3)

    for (key, matrix) in efizz
	sz = if matrix !== missing 
	    size(matrix)
	else
	    continue
	end
        if sz == (time_length, length(frequencies))  # Only process matrices with matching sizes
            band_means_matrix = Matrix{Float64}(undef, time_length, length(bands))

            for (band_idx, (band_name, (low, high))) in enumerate(bands)
                idx = findall(x -> x >= low && x <= high, frequencies)
                band_means_matrix[:, band_idx] = mean(matrix[:, idx], dims=2)
            end

            band_means_dict[key] = band_means_matrix
        else
            band_means_dict[key] = matrix  # keep other entries unchanged
        end
    end

    return band_means_dict
end


function convert_complex_to_real!(df::DataFrame)
    for col in names(df)
        if eltype(df[!, col]) <: Complex
            df[!, col] = real.(df[!, col])
        end
	# if eltype(df[!, col]) <: Float64
	#     df[!, col] = Float32.(df[!, col])
	# end
    end
    return df
end
function preprocess_dataframe!(df::DataFrame; zscore_all::Bool=false, quantile_normalize::Bool=false)
    for col in names(df)
        # If variance is zero, add small noise
        if var(df[!, col]) ≈ 0
            noise = randn(size(df[!, col])) * 1e-5
            df[!, col] .+= noise
	elseif zscore_all
	    if quantile_normalize
		lower, upper = quantile(df[!,col], [0.025, 0.975])
		df[!,col] = filter(x -> x>lower && x<upper, df[!,col])
	    end
            df[!, col] .= (df[!, col] .- mean(df[!, col])) ./ std(df[!, col])
        end
    end
    return df
end
function get_dataset(animal::String)
    filepath = "/Volumes/MATLAB-Drive/Shared/figures/midpattern=true/data/$(animal)_coherence_true.mat"
    efizzfile = "/Volumes/MATLAB-Drive/Shared/SingleDayExpt/$(animal)_direct/$(animal)spectralBehavior.mat";
    print("Reading matfile")
    mat = matread(filepath)
    println("...done")
    matefizz = matread(efizzfile)
    print("Loading data for animal $animal...")
    global U = mat["U"]
    global V = mat["V"]
    global R  = @. (U + V) / sqrt(2)
    global Ru = @. (U - V) / sqrt(2)
    global uv_time = mat["uv_time"]
    global spikeRateMatrix = mat["spikeRateMatrix"]
    efizz = matefizz["efizz"]
    global new_bands = generate_frequency_bands(vec(efizz["f"]), 3)
    efizz["phi_cos"] = cos.(efizz["phi"])
    efizz["phi_sin"] = sin.(efizz["phi"])
    efizza_var_to_use = ["S1", "S2", "Cavg", "wpli_avg", "phi_cos", "phi_sin"]
    global efizza_time = vec(efizz["t"])
    use_condensed = true 
    bands = generate_frequency_bands(vec(efizz["f"]), 3)
    global condensed_efizz = if use_condensed 
	compute_band_means_for_efizz(efizz)
    else
	efizz
    end
    time_indices = searchsortednearest.([efizza_time], uv_time)
    interpolated_efizz = Dict()
    println("Interpolating efizz data...")
    for field in efizza_var_to_use
        field_data = condensed_efizz[field]
        if size(field_data, 1) != length(efizza_time)
            println("Mismatch in dimensions for field: $field")
            continue
        end
	tmp  = field_data[time_indices, :]
	if typeof(tmp) <: Array{<:Complex}
	    tmp = real.(tmp)
	elseif typeof(tmp) <: Complex
	    tmp = real(tmp)
	end
	if size(tmp, 1) == 1
	    tmp = tmp[1,:,:,]
	end
        interpolated_data = tmp
        interpolated_efizz[field] = interpolated_data
    end
    println("Constructing dataframes ...")
    global df_R  = DataFrame(R[:,1:10], :auto)
    global df_Ru = DataFrame(Ru[:,1:10], :auto)
    rename!(df_R, [Symbol("R_", i) for i in 1:size(df_R, 2)])
    rename!(df_Ru, [Symbol("Ru_", i) for i in 1:size(df_Ru, 2)])
    global df_efizz = DataFrame()
    field = first(keys(interpolated_efizz))
    for field in efizza_var_to_use
        field_data = interpolated_efizz[field]
	freq_idx = 1
        for freq_idx in 1:size(field_data, 2)
	    if use_condensed
		column_name = Symbol("$(field)_$(collect(keys(bands))[freq_idx])")
	    else
		column_name = Symbol("$(field)_freq$(freq_idx)")
	    end
            tmp = field_data[:, freq_idx]
            if typeof(tmp) <: Complex || typeof(tmp) <: Vector{Complex}
                field_data[:, freq_idx] = real.(field_data[:, freq_idx])
            end
            if any(isnan, field_data[:, freq_idx])
                field_data[:, freq_idx] = fill(0.0, size(field_data[:, freq_idx]))
            end
            df_efizz[!, column_name] = field_data[:, freq_idx]
        end
    end
    df_efizz = convert_complex_to_real!(df_efizz)
    df_efizz = preprocess_dataframe!(df_efizz, zscore_all=true)
    return (;U, V, R, Ru, uv_time, spikeRateMatrix, df_R, df_Ru, df_efizz)
end
# function predict(test_data, coeffs)
#     # Assuming a linear combination for predictions
#     return test_data * coeffs'
# end
logistic(x) = 1. / (1. + exp(-x))
function get_test_indices(train_indices, total_data_size)
    return setdiff(1:total_data_size, train_indices)
end

# Load or initialize coefficients and models dictionaries
coef_dict = isfile(joinpath(folder, "efizza_coef_dict.jls")) ? deserialize(joinpath(folder, "efizza_coef_dict.jls")) : Dict{Tuple{String,String}, Any}()
models_dict = isfile(joinpath(folder, "efizza_models_dict.jls")) ? deserialize(joinpath(folder, "efizza_models_dict.jls")) : Dict{Tuple{String,String}, Any}()
inds_dict   = isfile(joinpath(folder, "efizza_inds_dict.jls")) ? 
    deserialize(joinpath(folder, "efizza_inds_dict.jls")) : Dict{String, Any}()
curr_animal = ""
stat_dict = isfile(joinpath(folder, "efizza_stat_dict.jls")) ? deserialize(joinpath(folder, "efizza_stat_dict.jls")) : Dict{Tuple{String,String,String}, Dict{String, Any}}()
predictions_df = isfile(joinpath(folder, "efizza_predictions_df.jls")) ? deserialize(joinpath(folder, "efizza_predictions_df.jls")) : DataFrame()
coeffs_df = isfile(joinpath(folder, "efizza_coeffs_df.jls")) ? deserialize(joinpath(folder, "efizza_coeffs_df.jls")) : DataFrame()
# if length(workers()) < n_chains
#     addprocs(n_chains - length(workers()))
#     @everywhere using TuringGLM, MCMCChains
# end
curranimal=""
animal = first(animals)

@showprogress "animal" for animal in animals
    if curranimal != animal
	global curranimal = animal
	global data = get_dataset(animal)
	print("Sampling for animal $animal ...")
    end
    @unpack U, V, R, Ru, uv_time, spikeRateMatrix, df_R, df_Ru, df_efizz = data

    df_all = DataFrames.hcat(df_R, df_Ru, df_efizz)
    if overwrite || !haskey(inds_dict, animal)
	if typeof(numsamprow) <: Int
	    inds = rand(1:size(df_all,1), numsamprow)
	elseif numsamprow == :numcols
	    inds = rand(1:size(df_all,1), size(df_all,2))
	else
	    @error "numsamprow must be an integer or :numcols"
	end
	inds_dict[animal] = inds
    else
	inds = inds_dict[animal]
    end
    # Separating into train and test data
    test_inds = get_test_indices(inds, size(df_all,1))
    train_data = df_all[inds, :]
    test_data = df_all[test_inds, :]
    train_data = DataFrame(convert.(Vector{Float32}, eachcol(train_data)), names(train_data))
    test_data = DataFrame(convert.(Vector{Float32}, eachcol(test_data)), names(test_data))

    # Train
    if usegpu
	# Convert columns to GPU arrays
	gpu_columns = OrderedDict()
	col = first(names(df_all))
	for col in names(df_all)
	    train_gpu_columns[col] = CuArray(df_all[!, col])
	    test_gpu_columns[col]  = CuArray(test_data[!, col])
	end
	# Convert dataframes to GPU arrays
	train_data = DataFrame(train_gpu_columns, names(train_data))
	test_data = DataFrame(test_gpu_columns, names(test_data))
    end
    predictors = names(df_efizz)
    println("Number of efizz predcitors: ", length(predictors))
    prog = Progress(length(vcat(names(df_R), names(df_Ru))), 1, "Processing columns...")
    println("Regressing each column...")
    target_cols = vcat(names(df_R)[1:2], names(df_Ru)[1:2])
    target_col = first(target_cols)
    @showprogress "targets" for target_col in target_cols

	if !haskey(coef_dict, (animal, target_col))
	    formula_str = "$(target_col) ~ " * join(predictors, " + ")
	    formula_expr = Meta.parse("formula = @formula($formula_str)")
	    formula = eval(formula_expr)
	    # md = turing_model(formula, train_data) # linear regression
            # let's use an error more appropriate for lfp
            md = turing_model(formula, train_data; model=TDist)
            # Custom Turing model with Gamma likelihood
            # Define a logpdf for SkewNormal distribution
            # function logpdf_skewnormal(x, location, scale, shape)
            #     z = (x - location) / scale
            #     return log(2) + logpdf(Normal(0, 1), z) + logcdf(Normal(0, 1), shape * z) - log(scale)
            # end
            #
            # @model function custom_model(y, X)
            #     # Priors for regression coefficients
            #     n, p = size(X)
            #     β ~ filldist(Normal(0, 2.5), p)
            #     
            #     # Priors for SkewNormal parameters
            #     location ~ Normal(0, 2)
            #     scale ~ Exponential(1)
            #     shape ~ Normal(0, 1)
            #     
            #     # Likelihood with Skewed Normal distribution
            #     for i in 1:n
            #         μ = dot(X[i, :], β)
            #         y[i] ~ logpdf_skewnormal(μ, location, scale, shape)
            #     end
            # end

	    tries = 5
	    while tries > -1
		try
		    # advi = ADVI(100, 1000)
		    # @time global q = vi(md, advi; optimizer=Flux.ADAM(0.1))
		    # @time global chns = sample(md, NUTS(), n_samples)
                    @time global chns = sample(md, NUTS(), n_samples)
                    # @time global chns = sample(custom_model(train_data[:, target_col], Matrix(train_data[:, predictors])), NUTS(), n_samples)
		    # global chns = @async sample(md, NUTS(), MCMCDistributed(),
		    # n_samples, n_chains; init_theta=init_vals)
		    chns = fetch(chns)
		    models_dict[(animal, target_col)] = chns
		    tries = -1
		    print("...done")
		catch e
		    println("Error encountered during sampling: ", e)
		    tries -= 1
		    continue
		end
	    end
	    coeffs = chns.value.data
	    println("Coefficients for target column $target_col: $coeffs")
	    coef_dict[(animal, target_col)] = coeffs
	end

	# TESTING
	paramset = "all"
        band_keys = (keys(bands)|>collect)|>x->replace.(string.(x), trades...)
	@showprogress "paramset" for paramset in ["all", "S1", "S2", "Cavg", "phi", band_keys..., ]
	    if haskey(stat_dict, (animal, target_col, paramset))
		continue
	    end
	    preds = if paramset == "all"
		predictors
	    else
		[x for x in predictors if occursin(paramset, x)]
	    end
	    if isempty(preds)
                continue
	    end
	    # Testing
	    stats = predict_and_compute_stats_v2(chns, test_data, preds, predictors, target_col)
	    # @unpack accuracy, r2, mse, predictions, prec, rec, f1_score = stats
	    println("Accuracy on test data: ", mean(stats.accuracy))
	    println("R2 on test data: ", mean(stats.r2))
	    println("MSE on test data: ", mean(stats.mse))
	    # Store the statistics into stat_dict
	    stat_dict[(animal, target_col, paramset)] = Dict(
		"Accuracy"=>stats.accuracy, "R2"=>stats.r2, "MSE" => stats.mse, 
		"Prec"=>stats.prec, "Rec"=>stats.rec, "F1_score"=>stats.f1_score
	    )
	    # Store the predictions into predictions_df
	    pred_df = DataFrame(
		animal = repeat([animal], length(test_inds)),
		commsub = repeat([target_col], length(test_inds)),
		time = uv_time[test_inds],
		actual = test_data[:, target_col],
		predicted = vec(mean(stats.predictions, dims=2)),
		precision = repeat([mean(stats.prec)], length(test_inds)),
		recall = repeat([mean(stats.rec)], length(test_inds)),
		f1_score = repeat([mean(stats.f1_score)], length(test_inds)),
		n_sample = repeat([n_samples], length(test_inds)),
		accuracy = repeat([mean(stats.accuracy)], length(test_inds)),
		r2 = repeat([mean(stats.r2)], length(test_inds)),
		mse = repeat([mean(stats.mse)], length(test_inds)),
		paramset = repeat([paramset], length(test_inds)),
	    )
	    append!(predictions_df, pred_df; cols=:union)
	end
    end

    if usegpu
	@everywhere gpu_columns = nothing
	@everywhere df_train = nothing
	@everywhere CUDA.reclaim()
	@everywhere GC.gc()
    end
    next!(prog)
    # Save the dictionaries to disk
    serialize(joinpath(folder, "efizza_coef_dict.jls"), coef_dict)
    serialize(joinpath(folder, "efizza_models_dict.jls"), models_dict)
    serialize(joinpath(folder, "efizza_inds_dict.jls"), inds_dict)
    serialize(joinpath(folder, "efizza_stat_dict.jls"), stat_dict)
    serialize(joinpath(folder, "efizza_predictions_df.jls"), predictions_df)
end



function replace_names(chns; mode=1)
    if mode  == 1
	existing_names = [nm for nm in names(chns) if occursin("β", string(nm))]
	new_names = propertynames(df_efizz)
    elseif mode == 2
	existing_names = propertynames(df_efizz)
	new_names = Symbol.("β" .* "_" .* string.(existing_names))
    elseif mode == 3
	existing_names = [nm for nm in names(chns) if occursin("β", string(nm))]
	new_names = Symbol.("β" .* "_" .* string.(propertynames(df_efizz)))
	new_names = Symbol.(replace.(string.(new_names), trades...))
    elseif mode == 4
	existing_names = Symbol.("β" .* "_" .* string.(propertynames(df_efizz)))
	# existing_names = Symbol.(replace.(string.(existing_names), r"phi_cos" => "phi-cos", r"phi_sin" => "phi-sin"))
	# existing_names = Symbol.(replace.(string.(existing_names), r"wpli_avg" => "wpli-avg"))
	new_names = Symbol.(replace.(string.(existing_names), trades...))
    else
	@error "Invalid mode"
    end
    replacements = OrderedDict{Symbol, Symbol}(existing_names[i] => new_names[i] for i in 1:length(existing_names))
    println("Replacing names...")
    println(replacements)
    replacenames(chns, replacements)
end

# Swap out each model's coefficients with the corresponding column names
for (key, chns) in models_dict
    models_dict[key] = replace_names(chns, mode=3)
end
serialize(joinpath(folder, "efizza_models_dict.jls"), models_dict)


stat_df = stat_dict_to_dataframe(stat_dict, true, :commsub)


animal = first(animals)
coefs_df = DataFrame()
if :efizz \notin propertynames(Main)
    data = get_dataset(animal)
    @unpack U, V, R, Ru, uv_time, spikeRateMatrix, df_R, df_Ru, df_efizz = data
end
println("Extracting coefficients for animal $animal...")
predictors = Symbol.(replace.(string.(propertynames(df_efizz)), trades...))
c_df = extract_coefficients(models_dict, string.(predictors); 
    target_prop=:commsub, splitter_props=[:Coefficient, :target_prop])
append!(coefs_df, c_df, cols=:union)
serialize(joinpath(folder, "efizza_coefs_df.jls"), coefs_df)

coefs_df.ValueAbs = abs.(coefs_df.Value)
coefs_df.orthostate = map(x->
begin
	if ismissing(x)
		missing
	elseif x == "R"
		"Aligned"
	elseif x=="Ru"
		"Orthogonal"
	end
end, coefs_df.target_prop1)
serialize(joinpath(folder, "efizza_coefs_df.jls"), coefs_df)

print("Adding lindist_bin_mid column...")
# Coefficients by animal and behavior
using PyCall
ENV["MPLBACKEND"] = "Agg"
@pyimport seaborn as sns
@pyimport pandas as pd
@pyimport matplotlib.pyplot as plt
@pyimport mplcyberpunk
plt.style.use("cyberpunk")
# Convert the DataFrame to a dictionary
Cdf = begin 
    Cdf = dropmissing(coefs_df, disallowmissing=true)
    Dict(names(Cdf) .=> eachcol(Cdf))
end
Cdf = pd.DataFrame(Cdf)
Cdf.to_csv(joinpath(folder, "efizza_coefs_df.csv"))

