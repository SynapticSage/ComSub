using Distributed: remotecall_eval
using Plots, StatsPlots
using Metrics
using Statistics
using DataFrames
using MAT
using Serialization
using Interpolations
import Turing
import MLJBase
using DataFrames
using ProgressMeter
using Distributed, Metrics
# import Threads
using SearchSortedNearest
using Infiltrator
using PyCall
import Turing, Flux
using Turing.Variational
@pyimport sklearn.metrics as skm
@pyimport warnings
warnings.filterwarnings("ignore")
@pyimport sklearn.preprocessing as skp
@everywhere using TuringGLM

function stat_dict_to_dataframe(data::Dict{Tuple{String, String, String}, Dict{String, Any}}, 
	explode=false, target_prop=:Behavior)
    # We will first figure out the unique stat keys from the first entry
    example_key = first(keys(data))
    stat_keys = keys(data[example_key])
    
    # Create columns for each of these stat keys
    columns = Dict{Symbol, Any}(:Animal => String[], target_prop => String[], :Type => String[])
    for stat_key in stat_keys
        columns[Symbol(stat_key)] = Vector{Float64}[]
    end

    all_rows = []
    # Extract data from the dictionary
    for ((animal, behavior, type), stats) in data
        if explode
            n = length(stats[first(stat_keys)])
            for i in 1:n
                row = Dict{Symbol, Any}(:Animal => animal, target_prop => behavior, :Type => type)
                for stat_key in stat_keys
                    row[Symbol(stat_key)] = stats[stat_key][i]
                end
                push!(all_rows, row)
            end
        else
            row = Dict{Symbol, Any}(:Animal => animal, target_prop => behavior, :Type => type)
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

# CONSTRUCT A DATAFRAME OF COEFFICIENTS
function extract_coefficients(models_dict::Dict, predictors::Vector{String};
	target_prop::Symbol=:Behavior, splitter_props=[:Coefficient])
	println("Current version")
    # Initialize an empty DataFrame
    coeffs_df = DataFrame(
        Animal = String[],
        target_prop = String[],
        Coefficient = String[],
        SampleNum = Int[],
        Value = Float64[]
    )
    # Iterate over the models_dict
    for ((animal, behavior), model) in models_dict
        # Extract α (intercept) samples
        α_samples = model[:α].data
        for (sample_num, α_val) in enumerate(α_samples)
            push!(coeffs_df, (animal, behavior, "Intercept", sample_num, α_val))
        end
        # Extract β (coefficients) samples
	@debug "Extracting coefficients for $animal $behavior"
	nms = filter(x->startswith(x,"β"), names(model).|>string)
	@assert(length(nms) == length(predictors))
	for (i,predictor) in enumerate(nms)
            # Extract the predictor's samples
            # predictor_name = "β[$predictor]"
	    predictor = string(predictor)
            predictor_name = predictor
	    @debug "...extracting $predictor_name"
            β_samples = model[predictor_name].data
            
            for (sample_num, β_val) in enumerate(β_samples)
                push!(coeffs_df, (animal, behavior, predictor, sample_num, β_val))
            end
        end
	@debug "...done!\n"
    end
    # Split the Coefficient column strings and create new columns
    if splitter_props == [:Coefficients]
	    coeffs_df[!, "Property1"] = Vector{Union{Missing,String}}(missing, size(coeffs_df, 1))
	    coeffs_df[!, "Property2"] = Vector{Union{Missing,String}}(missing, size(coeffs_df, 1))
	    for i in 1:size(coeffs_df, 1)
		# If the coefficient name is not "Intercept"
		if coeffs_df[i, :Coefficient] != "Intercept"
		    properties = split(coeffs_df[i, :Coefficient], "_")
		    coeffs_df[i, :Property1] = properties[1]
		    coeffs_df[i, :Property2] = properties[2]  # Extend this as needed for more properties
		else
		    coeffs_df[i, :Property1] = missing
		    coeffs_df[i, :Property2] = missing  # Extend this as needed for more properties
		end
	    end
    else
		for prop in splitter_props
			prop = string(prop)
			coeffs_df[!, prop * "1"] = Vector{Union{Missing,String}}(missing, size(coeffs_df, 1))
			coeffs_df[!, prop * "2"] = Vector{Union{Missing,String}}(missing, size(coeffs_df, 1))
			for i in 1:size(coeffs_df, 1)
				# If the coefficient name is not "Intercept"
				if occursin("_", coeffs_df[i, prop])
					properties = coeffs_df[i, prop]
					properties = replace.(properties, "β_" => "")
					properties = split(properties, "_")
					# @info "properties = $properties"
					coeffs_df[i, Symbol(prop * "1")] = properties[1]
					coeffs_df[i, Symbol(prop * "2")] = properties[2]  # Extend this as needed for more properties
				end
			end
		end
	end
    return coeffs_df
end

# For Bernoulli predictions, you might want to convert predictions to probabilities using logistic function
logistic(x) = 1. / (1. + exp(-x))

function get_test_indices(train_indices, total_data_size)
    return setdiff(1:total_data_size, train_indices)
end

animal = "JS21"
animals = ["JS21", "JS15", "JS14", "ZT2", "JS17", "JS13"]
