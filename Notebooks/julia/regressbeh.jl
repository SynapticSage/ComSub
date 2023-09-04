include("/Volumes/MATLAB-Drive/Shared/Notebooks/julia/regress_funcs.jl")
n_chains = 4;
# Number of samples
animals = ["JS21", "JS15", "JS14", "ZT2", "JS17", "JS13"]
n_samples = 500;
numsamprow = 50_000;
overwrite = true
Threads.nthreads() = n_chains;
# addprocs(n_chains-1)
folder = "/Volumes/MATLAB-Drive/Shared/figures/midpattern=true/data/"
global behavior_names = ["vel", "lindist", "accel", "trajbound", "leftright", "rewarded", "futureRewarded", "previousRewarded", "idphi", "future"]

if isfile(joinpath(folder, "coeffs_dict.jls"))
	coeffs_dict = deserialize(joinpath(folder, "coeffs_dict.jls"))
else
	coeffs_dict = Dict{Tuple{String,String}, Any}()
end
if isfile(joinpath(folder, "models_dict.jls"))
	models_dict = deserialize(joinpath(folder, "models_dict.jls"))
else
	models_dict = Dict{Tuple{String, String}, Any}()
end
if isfile(joinpath(folder, "inds_dict.jls"))
	inds_dict = deserialize(joinpath(folder, "inds_dict.jls"))
else
	inds_dict = Dict{Tuple{String, String}, Any}()
end
if isfile(joinpath(folder, "stat_dict.jls"))
	stat_dict = deserialize(joinpath(folder, "stat_dict.jls"))
else
	stat_dict = Dict{Tuple{String, String, String}, Dict{String, Any}}()
end
if isfile(joinpath(folder, "predictions_df.jls"))
	predictions_df = deserialize(joinpath(folder, "predictions_df.jls"))
else
	predictions_df = DataFrame()
end
if isfile(joinpath(folder, "coeffs_df.jls"))
	coefs_df = deserialize(joinpath(folder, "coeffs_df.jls"))
else
	coefs_df = DataFrame()
end

function get_dataset(animal)
	file = "/Volumes/MATLAB-Drive/Shared/figures/midpattern=true/data/$(animal)_coherence_true.mat";
	# Load data
	mat = matread(file);
	global U = mat["U"];
	global V = mat["V"];
	global R  = @. (U + V) / sqrt(2);
	global Ru = @. (U - V) / sqrt(2);
	global uv_time = mat["uv_time"];
	global spikeRateMatrix = mat["spikeRateMatrix"];
	global behavior_dict = mat["behavior_dict"];
	behavior_names = mat["behavior_names"];
	global behavior_time = vec(pop!(behavior_dict, "time"));
	pop!(behavior_dict, "X")
	pop!(behavior_dict, "Y")
	behavior_names = setdiff(behavior_names, ["time", "X", "Y"]);
	# For each behavior:
	time_indices = searchsortednearest.([behavior_time], vec(uv_time))
	for behavior in behavior_names
		if behavior ∉ keys(behavior_dict)
			continue
		end
		# Get the behavior data
		behavior_data = vec(behavior_dict[behavior])
		# Check if dimensions match
		if length(behavior_time) != length(behavior_data)
			println("Mismatch in dimensions for behavior: $behavior")
			println("Length of behavior_time: ", length(behavior_time))
			println("Length of behavior_data: ", length(behavior_data))
			continue
		end
		interpolated_behavior = behavior_data[time_indices]
		global behavior_dict[behavior] = interpolated_behavior
	end
	# Convert matrices U, V, and R to DataFrames for easier manipulation
	df_R  = DataFrame(R[:,1:10], :auto)
	df_Ru = DataFrame(Ru[:,1:10], :auto)
	# Rename columns for clarity
	rename!(df_R, [Symbol("R_", i) for i in 1:size(df_R, 2)])
	rename!(df_Ru, [Symbol("Ru_", i) for i in 1:size(df_Ru, 2)])
	(;U, V, R, Ru, uv_time, spikeRateMatrix, behavior_dict, behavior_names,
		df_R, df_Ru)
end
function get_df(pos...; behavior)
	""" Combine the position dataframes with the behavior dataframe """
	target_behavior = vec(behavior_dict[behavior])
	if any(isnan.(target_behavior))
		target_behavior = convert(Vector{Union{Missing, eltype(target_behavior)}}, target_behavior)
		target_behavior[isnan.(target_behavior)] .= missing
	end
	# Combine the dataframes
	global df = DataFrames.hcat(pos...)
	df[!,:Behavior].=target_behavior
	dropmissing!(df)
end

                              
#  . .     --.--          o     
# -+-+-      |  ,---.,---..,---.
# -+-+-      |  |    ,---|||   |
#  ` `       `  `    `---^``   '
                              

for animal in animals

	# Load data
	data = get_dataset(animal)
	@unpack U, V, R, Ru, uv_time, spikeRateMatrix, behavior_dict, behavior_names, df_R, df_Ru = data


	# For one behavior at a time:
	prog = Progress(length(behavior_names))
	@showprogress "behaviors" for behavior in behavior_names
		display(behavior)
		global bname = behavior
		if (animal,behavior) ∈ keys(coeffs_dict) && !overwrite
			next!(prog)
			continue
		end
		df = get_df(df_R, df_Ru; behavior=behavior)
		if isempty(df)
			next!(prog)
			continue
		end
		inds = rand(1:size(df,1), numsamprow)
		df = df[inds,:]
		if length(unique(df.Behavior)) == 1
			next!(prog)
			continue
		end
		# Fit GLM
		predictors = Symbol.(vcat(names(df_R), names(df_Ru)))
		formula_str = "Behavior ~ " * join(predictors, " + ")
		# Convert the string to an expression
		formula_expr = Meta.parse("formula = @formula($formula_str)")
		target_behavior = df[:, :Behavior]
		ubeh = unique(target_behavior)
		if length(ubeh) == 2
			model = Bernoulli;
			init_vals = []
		elseif length(ubeh) > 2 && length(ubeh) < 20
			model = NegativeBinomial;
			init_vals = []
		else
			model = TDist;
			init_vals = [(:ν, 1.0), (:σ, 4.0)]  # or any other positive value
		end
		formula = eval(formula_expr)
		md      = turing_model(formula, df; model)
		# Evaluate the expression
		# Sample using the NUTS sampler and parallel chains
		# @time global chns = sample(md, NUTS(), Turing.MCMCThreads(), n_samples, n_chains)

		tries = 5
		while tries > -1
			try
				@time global chns = sample(md, NUTS(), n_samples; init_theta=init_vals)
			catch
				tries -= 1
				continue
			end
			models_dict[animal,behavior] = chns
			tries = -1
		end
		# Extract and interpret results
		# coeffs = coef(chns)
		coeffs = chns.value.data;
		println("Coefficients for behavior $behavior: $coeffs")
			
		# Store coefficients and models in the dictionaries
		coeffs_dict[animal,behavior] = coeffs
		inds_dict[animal,behavior]   = inds
		next!(prog)
	end

	# Save the dictionaries to disk
	serialize(joinpath(folder, "coeffs_dict.jls"), coeffs_dict)
	serialize(joinpath(folder, "models_dict.jls"), models_dict)
	serialize(joinpath(folder, "inds_dict.jls"), inds_dict)
end

                             
#  . .     --.--          |    
# -+-+-      |  ,---.,---.|--- 
# -+-+-      |  |---'`---.|    
#  ` `       `  `---'`---'`---'
                             

"""
	predict_and_compute_stats(model::Chains, 
		test_data::DataFrame, 
		paramrange::Union{Vector,UnitRange}=nothing)
Predicts on test data and computes summary statistics.
- `model` is the chains object returned by Turing
- `test_data` is a DataFrame containing the test data
- `paramrange` is a vector or range of indices of the parameters to use for prediction
"""
function predict_and_compute_stats_v2(model,
	test_data, predictors;
	paramrange=nothing)
	print("Predicting on test data...")
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
    @assert(size(X_test, 2) == size(β, 2))
    @assert(size(X_test, 1) == size(test_data, 1))
	predictions = zeros(size(X_test,1), size(α,1))
	for i in eachindex(α)
		predictions[:,i] = α[i] .+ (X_test * β[i, :].data)
	end
    # Compute summary stats
    true_values = test_data[:, :Behavior]  # Assuming "Behavior" is your target column
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

prog = Progress(length(Iterators.product(animals, behavior_names)))
animal = first(animals)
@showprogress "animal" for animal in animals
	data = get_dataset(animal)
	@unpack U, V, R, Ru, uv_time, spikeRateMatrix, behavior_dict, behavior_names, df_R, df_Ru = data
	behavior = first(behavior_names)
	@showprogress "testing behaviors" for behavior in behavior_names
		# Load data
		if (animal,behavior) ∉ keys(coeffs_dict) || 
			(animal,behavior) ∉ keys(models_dict)  ||
			(animal,behavior) ∉ keys(inds_dict)
			@warn "Skipping $(animal) $(behavior)"
			if (animal,behavior) ∉ keys(models_dict)
				@warn "Skipping $(animal) $(behavior)"
			end
			next!(prog)
			continue
		end
		@info "engaging with $(animal) $(behavior)"
		# sleep(1)
		@infiltrate
		inds   = inds_dict[animal,behavior]
		model     = models_dict[animal,behavior]
		predictors = Symbol.(vcat(names(df_R), names(df_Ru)))
		df = get_df(df_R, df_Ru; behavior=behavior)
		# Identify testing set
		test_indices = get_test_indices(inds, size(df, 1))
		# Filter the dataframe to only include the test set
		df_test = df[test_indices, :]
		# Predict and compute stats
		@showprogress "paramset" for paramset in ["R", "Ru", "all"]
			if (animal,behavior,paramset) ∈ keys(stat_dict)
				next!(prog)
				println("Skipping $(animal) $(behavior) $(paramset)")
				continue
			end
			@info "engaging with $(animal) $(behavior) $(paramset)"
			propnames = if paramset == "R"
				propertynames(df_R)
			elseif paramset == "Ru"
				propertynames(df_Ru)
			else
				predictors
			end
			@time dat = predict_and_compute_stats_v2(model, df_test, propnames; paramrange= 1:length(propnames))
			@unpack accuracy, r2, mse, predictions, prec, rec, f1_score = dat
			predictions = vec(mean(predictions, dims=2))
			# Store the statistics into stat_dict
			stat_dict[(animal, behavior, paramset)] = Dict(
				"Accuracy"=>dat.accuracy, 
				"R2"=>dat.r2,
				"MSE" => dat.mse,
				"Prec"=>dat.prec,
				"Rec"=>dat.rec,
				"F1_score"=>dat.f1_score)
			# Store the predictions into predictions_df
			pred_df = DataFrame(
				animal = repeat([animal], length(test_indices)),
				behavior = repeat([behavior], length(test_indices)),
				time = uv_time[test_indices],
				actual = df_test[:, :Behavior],
				predicted = vec(mean(dat.predictions, dims=2)),
				n_sample = repeat([length(test_indices)], length(test_indices)),
				accuracy = repeat([mean(dat.accuracy)], length(test_indices)),
				r2 = repeat([mean(dat.r2)], length(test_indices)),
				mse = repeat([mean(dat.mse)], length(test_indices)),
				precision = repeat([mean(dat.prec)], length(test_indices)),
				recall = repeat([mean(dat.rec)], length(test_indices)),
				f1_score = repeat([mean(dat.f1_score)], length(test_indices)),
				paramset = repeat(["all"], length(test_indices))
			)
			append!(predictions_df, pred_df)
		end
	end
	serialize(joinpath(folder, "stat_dict.jls"), stat_dict)
	serialize(joinpath(folder, "predictions_df.jls"), predictions_df)
end


# Usage
stat_df = stat_dict_to_dataframe(stat_dict, true)  # Set true to explode the vectors

                                                          
#  . .     ,---.          ,---.         |         |         
# -+-+-    |    ,---.,---.|__.     ,---.|    ,---.|--- ,---.
# -+-+-    |    |   ||---'|        |   ||    |   ||    `---.
#  ` `     `---'`---'`---'`        |---'`---'`---'`---'`---'
#                                  |                        


animal = first(animals)
coefs_df = DataFrame()
@showprogress "animal coef dict" for animal in animals
	data = get_dataset(animal)
	@unpack U, V, R, Ru, uv_time, spikeRateMatrix, behavior_dict, behavior_names, df_R, df_Ru = data
	predictors = Symbol.(vcat(names(df_R), names(df_Ru)))
	c_df = extract_coefficients(models_dict, string.(predictors))
	append!(coefs_df, c_df, cols=:union)
	serialize(joinpath(folder, "coefs_df.jls"), coefs_df)
end
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
end, coefs_df.Property1)
serialize(joinpath(folder, "coefs_df.jls"), coefs_df)

