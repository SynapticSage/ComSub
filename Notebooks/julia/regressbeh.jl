using Plots, StatsPlots
using Metrics
using Statistics
using DataFrames
using MAT
using Serialization
using Interpolations
import Turing
using DataFrames
using ProgressMeter
using Distributed, Metrics
# import Threads
using SearchSortedNearest
using Infiltrator
@everywhere using TuringGLM
n_chains = 4;
# Number of samples
n_samples = 500;
numsamprow = 50_000;
Threads.nthreads() = n_chains;
# addprocs(n_chains-1)
folder = "/Volumes/MATLAB-Drive/Shared/figures/midpattern=true/data/"
plotfolder=joinpath(splitpath(folder)[1:end-1]..., "julia_regressbeh")
if !isdir(plotfolder)
	mkdir(plotfolder)
end

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
	predictions_df = DataFrame(animal = String[], behavior = String[], time = Float64[], actual = Float64[], predicted = Float64[], n_sample = Int[],
		accuracy = Float64[], r2 = Float64[], mse = Float64[], paramset = String[])
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
	global behavior_names = ["vel", "lindist", "accel", "trajbound", "leftright", "rewarded", "futureRewarded", "previousRewarded", "idphi", "future"]
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

animals = ["JS21", "JS15", "JS14", "ZT2", "JS17", "JS13"]

for animal in animals

	# Load data
	data = get_dataset(animal)
	@unpack U, V, R, Ru, uv_time, spikeRateMatrix, behavior_dict, behavior_names, df_R, df_Ru = data


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

	# For one behavior at a time:
	prog = Progress(length(behavior_names))
	@showprogress "behaviors" for behavior in behavior_names
		display(behavior)
		global bname = behavior
		if (animal,behavior) ∈ keys(coeffs_dict)
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
end



# Function to predict based on test data and coefficients
function predict(test_data, coeffs)
    # Assuming a linear combination for predictions
    return test_data * coeffs'
end

# For Bernoulli predictions, you might want to convert predictions to probabilities using logistic function
logistic(x) = 1. / (1. + exp(-x))

function get_test_indices(train_indices, total_data_size)
    return setdiff(1:total_data_size, train_indices)
end

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
	using Metrics
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
	binary_preds = round.(predictions)
	accuracy = vec(mean(binary_preds .== true_values, dims=1))
	precision = vec(precision(binary_preds, true_values))
	recall = vec(recall(binary_preds, true_values))
	f1_score = 2 .* ((precision .* recall) ./ (precision .+ recall))
	# MSE
	mse = vec(mean(residuals .^ 2, dims=1))
	print("Accuracy: ", mean(accuracy))
	print("R2: ", mean(r2))
	print("MSE: ", mean(mse))
    return (;accuracy, r2, mse, predictions, precision, recall, f1_score)
end

stat_dict = Dict{Tuple{String, String, String}, Dict{String, Any}}()
predictions_df = DataFrame(animal = String[], behavior = String[], time = Float64[], actual = Float64[], predicted = Float64[], n_sample = Int[],
	accuracy = Float64[], r2 = Float64[], mse = Float64[], paramset = String[])
prog = Progress(length(Iterators.product(animals, behavior_names)))
@showprogress "animal" for animal in animals
	data = get_dataset(animal)
	@unpack U, V, R, Ru, uv_time, spikeRateMatrix, behavior_dict, behavior_names, df_R, df_Ru = data
	@showprogress "testing behaviors" for behavior in behavior_names
		display(behavior)
		# Load data
		if (animal,behavior) ∈ keys(stat_dict) || (animal,behavior) ∉ keys(coeffs_dict)
			next!(prog)
			continue
		end
		inds   = inds_dict[animal,behavior]
		model     = models_dict[animal,behavior]
		predictors = Symbol.(vcat(names(df_R), names(df_Ru)))
		df = get_df(df_R, df_Ru; behavior=behavior)
		# Identify testing set
		test_indices = get_test_indices(inds, size(df, 1))
		# Filter the dataframe to only include the test set
		df_test = df[test_indices, :]
		# Predict and compute stats
		begin
			@time accuracy, r2, mse, predictions = predict_and_compute_stats_v2(model, df_test, predictors)
			predictions = vec(mean(predictions, dims=2))
			# Store the statistics into stat_dict
			stat_dict[(animal, behavior, "all")] = Dict(
				"Accuracy"=>accuracy, 
				"R2"=>r2,
				"MSE" => mse)
			# Store the predictions into predictions_df
			pred_df = DataFrame(
				animal = repeat([animal], length(test_indices)),
				behavior = repeat([behavior], length(test_indices)),
				time = uv_time[test_indices],
				actual = df_test[:, :Behavior],
				predicted = predictions,
				n_sample = repeat([length(test_indices)], length(test_indices)),
				accuracy = repeat([mean(accuracy)], length(test_indices)),
			r2 = repeat([mean(r2)], length(test_indices)),
				mse = repeat([mean(mse)], length(test_indices)),
				paramset = repeat(["all"], length(test_indices))
			)
			append!(predictions_df, pred_df)
		end
		# For R variables
		begin
			@time accuracy, r2, mse, predictions = predict_and_compute_stats_v2(model, df_test, propertynames(df_R);paramrange= 1:size(df_R,2))
			predictions = vec(mean(predictions, dims=2))
			# Store the statistics into stat_dict
			stat_dict[(animal, behavior, "R")] = Dict(
				"Accuracy"=>accuracy, 
				"R2"=>r2,
				"MSE" => mse)
			# Store the predictions into predictions_df
			pred_df = DataFrame(
				animal = repeat([animal], length(test_indices)),
				behavior = repeat([behavior], length(test_indices)),
				time = uv_time[test_indices],
				actual = df_test[:, :Behavior],
				predicted = predictions,
				n_sample = repeat([length(test_indices)], length(test_indices)),
				accuracy = repeat([mean(accuracy)], length(test_indices)),
			r2 = repeat([mean(r2)], length(test_indices)),
			mse = repeat([mean(mse)], length(test_indices)),
				paramset = repeat(["R"], length(test_indices))
			)
			append!(predictions_df, pred_df)
		end
		# For Ru variables
		begin
			@time accuracy, r2, mse, predictions = predict_and_compute_stats_v2(model, df_test, propertynames(df_Ru);paramrange= 1:size(df_Ru,2))
			predictions = vec(mean(predictions, dims=2))
			# Store the statistics into stat_dict
			stat_dict[(animal, behavior, "Ru")] = Dict(
				"Accuracy"=>accuracy, 
				"R2"=>r2,
				"MSE" => mse)
			# Store the predictions into predictions_df
			pred_df = DataFrame(
				animal = repeat([animal], length(test_indices)),
				behavior = repeat([behavior], length(test_indices)),
				time = uv_time[test_indices],
				actual = df_test[:, :Behavior],
				predicted = predictions,
				n_sample = repeat([length(test_indices)], length(test_indices)),
				accuracy = repeat([mean(accuracy)], length(test_indices)),
			r2 = repeat([mean(r2)], length(test_indices)),
			mse = repeat([mean(mse)], length(test_indices)),
				paramset = repeat(["Ru"], length(test_indices))
			)
			append!(predictions_df, pred_df)
		end
		# Progress
		next!(prog)
	end
	serialize(joinpath(folder, "stat_dict.jls"), stat_dict)
	serialize(joinpath(folder, "predictions_df.jls"), predictions_df)
end


"""
	plot_predictions(animal, behavior, paramset=nothing)
Plots the predictions for a given animal and behavior.
- `animal` is the animal name
- `behavior` is the behavior name
- `paramset` is a vector or range of indices of the parameters to use for prediction
Output: a plot object
"""
function plot_predictions(animal, behavior, paramset=nothing)
	mdl = models_dict[animal, behavior]
	if paramset === nothing
		nms = [n for n in names(mdl) if contains(string(n), "β") || contains(string(n), "α")]
		paramset = 1:length(nms)
	end
	p=plot(mdl[paramset], legend=false)
	# Assuming you have 33 rows and 2 columns
	for i in paramset
		subplot = p[i, 2]  # get the subplot in the 2nd column for each row
		subplot2 = p[i, 1]  # get the subplot in the 1st column for each row
		q = subplot.series_list[1]  # get the plot object
		q2 = subplot2.series_list[1]  # get the plot object
		if i > 1 && i <= 11
			q.plotattributes[:fillcolor] = :blue  # set the color attribute of the plot
			q2.plotattributes[:fillcolor] = :blue  # set the color attribute of the plot
		elseif i > 11 && i <= 22
			q.plotattributes[:fillcolor] = :red
			q2.plotattributes[:fillcolor] = :red
		else
			q.plotattributes[:fillcolor] = :blue
			q2.plotattributes[:fillcolor] = :blue
		end
		vline!(subplot, [0], color=:red)  # add a vertical line at x=0
	end
	plot!(p, link=:x)
	p
end



function stat_dict_to_dataframe(data::Dict{Tuple{String, String, String}, Dict{String, Any}}, explode=false)
    # We will first figure out the unique stat keys from the first entry
    example_key = first(keys(data))
    stat_keys = keys(data[example_key])
    
    # Create columns for each of these stat keys
    columns = Dict{Symbol, Any}(:Animal => String[], :Behavior => String[], :Type => String[])
    for stat_key in stat_keys
        columns[Symbol(stat_key)] = Vector{Float64}[]
    end

    all_rows = []

    # Extract data from the dictionary
    for ((animal, behavior, type), stats) in data
        if explode
            n = length(stats[first(stat_keys)])
            for i in 1:n
                row = Dict{Symbol, Any}(:Animal => animal, :Behavior => behavior, :Type => type)
                for stat_key in stat_keys
                    row[Symbol(stat_key)] = stats[stat_key][i]
                end
                push!(all_rows, row)
            end
        else
            row = Dict{Symbol, Any}(:Animal => animal, :Behavior => behavior, :Type => type)
            for stat_key in stat_keys
                row[Symbol(stat_key)] = stats[stat_key]
            end
            push!(all_rows, row)
        end
    end

    # Convert rows to DataFrame
    df = DataFrame(all_rows)
    return df
end

# Usage
stat_df = stat_dict_to_dataframe(stat_dict, true)  # Set true to explode the vectors

function generate_heatmap(df::DataFrame, metric::Symbol; type=nothing, diff_type=nothing, kwargs...)
    
    # If diff_type is provided, ensure type is also provided
    if diff_type !== nothing && type === nothing
        throw(ArgumentError("If diff_type is provided, type must also be specified"))
    end
    
    # Create matrix for the heatmap
    matrix_data = Matrix{Float64}(undef, length(unique(df.Animal)), length(unique(df.Behavior)))
    
    for (i, animal) in enumerate(unique(df.Animal))
        for (j, behavior) in enumerate(unique(df.Behavior))
            
            # If diff_type is not provided, use the original type behavior
            if diff_type === nothing
                subset = type === nothing ? df[(df.Animal .== animal) .& (df.Behavior .== behavior), :] : 
                                             filter(row -> row.Type == type, df[(df.Animal .== animal) .& (df.Behavior .== behavior), :])
                matrix_data[i, j] = mean(subset[:, metric])
                
            # If diff_type is provided, compute the difference
            else
                subset_type = filter(row -> row.Type == type, df[(df.Animal .== animal) .& (df.Behavior .== behavior), :])
                subset_diff_type = filter(row -> row.Type == diff_type, df[(df.Animal .== animal) .& (df.Behavior .== behavior), :])
                
                matrix_data[i, j] = mean(subset_type[:, metric]) - mean(subset_diff_type[:, metric])
            end
        end
    end
    
    title_str = diff_type === nothing ? "Heatmap of Mean $metric" : "Heatmap of Difference in $metric between $type and $diff_type"
    
    defaults = (;xlabel="Animal", ylabel="Behavior", title=title_str, color=:viridis, aspect_ratio=:auto)
    kw = merge(defaults, kwargs)
    
    # Generate heatmap
    heatmap(unique(df.Animal), unique(df.Behavior), matrix_data'; kw...)
end
generate_heatmap(stat_df, :Accuracy, color=:balance)
savefig(joinpath(plotfolder, "heatmap_Accuracy.png"))
savefig(joinpath(plotfolder, "heatmap_Accuracy.pdf"))
generate_heatmap(stat_df, :R2, color=:balance, clim=(-1,1))
savefig(joinpath(plotfolder, "heatmap_R2.png"))
savefig(joinpath(plotfolder, "heatmap_R2.pdf"))
sort!(stat_df, [:Type,:Behavior,:Animal])
plot(generate_heatmap(stat_df, :Accuracy, type="R", color=:balance, title="CommSub"),
	generate_heatmap(stat_df, :Accuracy, type="Ru", color=:balance, title="Ortthogonal Commsub"), size=(1000,500))
savefig(joinpath(plotfolder, "heatmap_Accuracy_R_Ru.png"))
savefig(joinpath(plotfolder, "heatmap_Accuracy_R_Ru.pdf"))
plot(generate_heatmap(stat_df, :R2, type="R", color=:balance, title="CommSub"),
	generate_heatmap(stat_df, :R2, type="Ru", color=:balance, title="Orthogonal Commsub"), size=(1000,500), clim=(-0.5, 0.5))
savefig(joinpath(plotfolder, "heatmap_R2_R_Ru.png"))
savefig(joinpath(plotfolder, "heatmap_R2_R_Ru.pdf"))
plot(
generate_heatmap(stat_df, :Accuracy, type="R", diff_type="Ru", title="ComSub-Orthogonal", color=:Reds),
	generate_heatmap(stat_df, :Accuracy, type="all", diff_type="R", title="All-R: Orthog Improve Commsub", color=:Purples, clim=(0,0.125)),
size=(1000,500), 
)
savefig(joinpath(plotfolder, "heatmap_Accuracy_R_Ru_diff.png"))
savefig(joinpath(plotfolder, "heatmap_Accuracy_R_Ru_diff.pdf"))

summ = combine(groupby(stat_df, [:Behavior, :Type]),
	:Accuracy => mean,
	:R2 => mean,
	:MSE => mean,
	:Accuracy => (x -> quantile(x, 0.025)) => :Accuracy_lower,
	:Accuracy => (x -> quantile(x, 0.975)) => :Accuracy_upper,
	:R2 => (x -> (quantile(x, 0.025))) => :R2_lower,
	:R2 => (x -> quantile(x, 0.975)) => :R2_upper,
	:MSE => (x -> quantile(x, 0.025)) => :MSE_lower,
	:MSE => (x -> quantile(x, 0.975)) => :MSE_upper,
	
)

function generate_barplot(df::DataFrame)
    # Create a grouped bar plot for each metric
    p1 = groupedbar(df.Behavior, df[:, :Accuracy_mean], group=df.Type, 
                    yerr=(df.Accuracy_mean .- df.Accuracy_lower, df.Accuracy_upper .- df.Accuracy_mean), 
                    title="Accuracy", legend=:topright)
	hline!([0], color=:black, linestyle=:dash)
    p2 = groupedbar(df.Behavior, df[:, :R2_mean], group=df.Type, 
                    yerr=(df.R2_mean .- df.R2_lower, df.R2_upper.-df.R2_mean), 
		title="R2", legend=:topright, ylim=(-0.5,0.75))
	hline!([0], color=:black, linestyle=:dash)
    # Combine the subplots
	plot(p1, p2, layout=(2, 1), size=0.8.*(1000,800), bbox=(0.5, 0.8))
end
generate_barplot(summ)
savefig(joinpath(plotfolder, "barplot_stats_by_behavior_type.png"))
savefig(joinpath(plotfolder, "barplot_stats_by_behavior_type.pdf"))

summ = combine(groupby(stat_df, [:Behavior, :Type, :Animal]),
	:Accuracy => mean,
	:R2 => mean,
	:MSE => mean,
	:Accuracy => (x -> quantile(x, 0.025)) => :Accuracy_lower,
	:Accuracy => (x -> quantile(x, 0.975)) => :Accuracy_upper,
	:R2 => (x -> (quantile(x, 0.025))) => :R2_lower,
	:R2 => (x -> quantile(x, 0.975)) => :R2_upper,
	:MSE => (x -> quantile(x, 0.025)) => :MSE_lower,
	:MSE => (x -> quantile(x, 0.975)) => :MSE_upper,
	
)


function generate_animal_barplots(df::DataFrame)
    # Get unique animals
    animals = unique(df.Animal)
    plots = []
    for animal in animals
        # Filter dataframe by animal
        df_animal = filter(row -> row.Animal == animal, df)
        
        # Create subplots for the current animal
        p1 = groupedbar(df_animal.Behavior, df_animal[:, :Accuracy_mean], group=df_animal.Type,
                        yerr=(df_animal.Accuracy_mean .- df_animal.Accuracy_lower, df_animal.Accuracy_upper .- df_animal.Accuracy_mean),
                        title="Accuracy", legend=false,
						ylims=(0,1))
        
        p2 = groupedbar(df_animal.Behavior, df_animal[:, :R2_mean], group=df_animal.Type,
                        yerr=(df_animal.R2_mean .- df_animal.R2_lower, df_animal.R2_upper .- df_animal.R2_mean),
                        title="R2", legend=false, ylim=(-0.1,0.4))
        
        
        # Combine the subplots for the current animal
        animal_plot = plot(p1, p2, layout=(1, 2), title=animal)
        push!(plots, animal_plot)
    end
    
    # Combine all the plots
	final_plot = plot(plots..., layout=(length(animals), 1), size=1.1.*(1100, 150 * length(animals)), bbox=(0.2, 0.99), xrotation=15, fontsize=14)
    return final_plot
end
generate_animal_barplots(summ)
savefig(joinpath(plotfolder, "barplot_stats_by_behavior_type_animal.png"))
savefig(joinpath(plotfolder, "barplot_stats_by_behavior_type_animal.pdf"))

plots = Dict()
group = first(groupby(stat_df, [:Behavior, :Type, :Animal]))
for group in groupby(stat_df,  [:Behavior, :Type, :Animal])
	if isempty(group)
		continue
	end
	try
		p=plot(
			histogram(group.Accuracy, title="$(group.Animal[1]) $(group.Behavior[1]) $(group.Type[1])", xlabel="Accuracy", bins=25),
			histogram(group.R2, title="$(group.Animal[1]) $(group.Behavior[1]) $(group.Type[1])", xlabel="R2"),
		)
		push!(plots, (group.Behavior[1], group.Type[1], group.Animal[1]) => p)
		savefig(joinpath(plotfolder, "stat_$(group.Animal[1])_$(group.Behavior[1])_$(group.Type[1]).png"))
		savefig(joinpath(plotfolder, "stat_$(group.Animal[1])_$(group.Behavior[1])_$(group.Type[1]).pdf"))
	catch ME
	end
end


