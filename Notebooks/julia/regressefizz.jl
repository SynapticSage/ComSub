using DataFrames, Statistics, LinearAlgebra, Serialization
using Interpolations
using MAT
using Plots
unicodeplots()
using Distributed
using ProgressMeter
using DataStructures: OrderedDict

import Turing
using Turing.Variational
using SearchSortedNearest

animal = "JS21"
folder = "/Volumes/MATLAB-Drive/Shared/figures/midpattern=true/data/"
file = "/Volumes/MATLAB-Drive/Shared/figures/midpattern=true/data/$(animal)_coherence_true.mat";


usegpu = true
n_chains = 16;
# Number of samples
n_samples = Int(round(150/n_chains));
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


function compute_band_means(efizz_data, frequencies)
    bands = Dict(
        "delta" => (0.5, 4.0),
        "theta" => (6.0, 12.0),
        "beta" => (13.0, 30.0),
        "low_gamma" => (30.0, 50.0),
        "mid_gamma" => (50.0, 80.0),
        "high_gamma" => (80.0, 100.0),
        "epsilon" => (100.0, 120.0),
        "ripple" => (120.0, 160.0)
    )
    band_means = Dict()
    for (band_name, (low, high)) in bands
        idx = findall(x -> x >= low && x <= high, frequencies)
        band_means[band_name] = mean(efizz_data[:, idx], dims=2)
    end
    return band_means
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
function preprocess_dataframe!(df::DataFrame; zscore_all::Bool=false)
    for col in names(df)
        # If variance is zero, add small noise
        if var(df[!, col]) ≈ 0
            noise = randn(size(df[!, col])) * 1e-5
            df[!, col] .+= noise
	elseif zscore_all
            df[!, col] .= (df[!, col] .- mean(df[!, col])) ./ std(df[!, col])
        end
    end
    return df
end
function get_dataset(animal::String)
    filepath = "/Volumes/MATLAB-Drive/Shared/figures/midpattern=true/data/$(animal)_coherence_true.mat"
    efizzfile = "/Volumes/MATLAB-Drive/Shared/SingleDayExpt/$(animal)_direct/$(animal)spectralBehavior.mat";
    mat = matread(filepath)
    matefizz = matread(efizzfile)
    global U = mat["U"]
    global V = mat["V"]
    global R  = @. (U + V) / sqrt(2)
    global Ru = @. (U - V) / sqrt(2)
    global uv_time = mat["uv_time"]
    global spikeRateMatrix = mat["spikeRateMatrix"]
    global efizz = matefizz["efizz"]
    efizz["phi_cos"] = cos.(efizz["phi"])
    efizz["phi_sin"] = sin.(efizz["phi"])
    efizz_var_to_use = ["S1", "S2", "Cavg", "wpli_avg", "phi_cos", "phi_sin"]
    global efizz_time = vec(efizz["t"])
    time_indices = searchsortednearest.([efizz_time], uv_time)
    interpolated_efizz = Dict()
    println("Interpolating efizz data...")
    for field in efizz_var_to_use
        field_data = efizz[field]
        if size(field_data, 1) != length(efizz_time)
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
    println("Constructing dataframes...")
    global df_R  = DataFrame(R[:,1:10], :auto)
    global df_Ru = DataFrame(Ru[:,1:10], :auto)
    rename!(df_R, [Symbol("R_", i) for i in 1:size(df_R, 2)])
    rename!(df_Ru, [Symbol("Ru_", i) for i in 1:size(df_Ru, 2)])
    global df_efizz = DataFrame()
    for field in efizz_var_to_use
        field_data = interpolated_efizz[field]
        for freq_idx in 1:size(field_data, 2)
            column_name = Symbol("$(field)_freq$(freq_idx)")
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


# Load or initialize coefficients and models dictionaries
coeffs_dict = isfile(joinpath(folder, "efizz_coeffs_dict.jls")) ? deserialize(joinpath(folder, "coeffs_dict.jls")) : Dict{String, Any}()
models_dict = isfile(joinpath(folder, "efizz_models_dict.jls")) ? deserialize(joinpath(folder, "models_dict.jls")) : Dict{String, Any}()
inds_dict = isfile(joinpath(folder, "efizz_inds_dict.jls")) ? deserialize(joinpath(folder, "inds_dict.jls")) : Dict{String, Any}()
curr_animal = ""

if length(workers()) < n_chains
    addprocs(n_chains - length(workers()))
    @everywhere using TuringGLM, MCMCChains
end
animals = ["JS21", "JS15", "JS14", "ZT2", "JS17", "JS13"]
animal = first(animals)
@showprogress "animal" for animal in animals

    global curranimal = animal
    global data = get_dataset(animal)
    @unpack U, V, R, Ru, uv_time, spikeRateMatrix, df_R, df_Ru, df_efizz = data
    print("Sampling for animal $animal...")

    df_all = DataFrames.hcat(df_R, df_Ru, df_efizz)
    if typeof(numsamprow) <: Int
	inds = rand(1:size(df_all,1), numsamprow)
    elseif numsamprow == :numcols
	inds = rand(1:size(df_all,1), size(df_all,2))
    else
	@error "numsamprow must be an integer or :numcols"
    end
    df_all = df_all[inds, :]
    df_all = DataFrame(convert.(Vector{Float32}, eachcol(df_all)), names(df_all));
    if usegpu
	# Convert columns to GPU arrays
	gpu_columns = OrderedDict()
	col = first(names(df_all))
	for col in names(df_all)
	    gpu_columns[col] = CuArray(df_all[!, col])
	end
	df_all = DataFrame(gpu_columns)
    end

    predictors = names(df_efizz)
    prog = Progress(length(vcat(names(df_R), names(df_Ru))), 1, "Processing columns...")

    println("Regressing each column...")
    target_col = first(vcat(names(df_R), names(df_Ru)))
    @showprogress "targets" for target_col in vcat(names(df_R), names(df_Ru))

        formula_str = "$(target_col) ~ " * join(predictors, " + ")
        formula_expr = Meta.parse("formula = @formula($formula_str)")
        formula = eval(formula_expr)
	model = TDist
	init_vals = [(:ν, 1.0), (:σ, 4.0)]
	md = turing_model(formula, df_all; model=model)

        tries = 5
        while tries > -1
            try
		# advi = ADVI(1_000, 1000)
		# @time global q = vi(md, advi; optimizer=Flux.ADAM(0.1))
                @time global chns = sample(md, NUTS(), n_samples; init_theta=init_vals)
		# global chns = @async sample(md, NUTS(), MCMCDistributed(),
		    # n_samples, n_chains; init_theta=init_vals)
		chns = fetch(chns)
                models_dict[(animal, target_col)] = chns
                tries = -1
            catch e
                println("Error encountered during sampling: ", e)
                tries -= 1
                continue
            end
        end

        coeffs = chns.value.data
        println("Coefficients for target column $target_col: $coeffs")
        coeffs_dict[(animal, target_col)] = coeffs
        inds_dict[(animal, target_col)] = inds
	if usegpu
	    @everywhere gpu_columns = nothing
	    @everywhere df_all = nothing
	    @everywhere CUDA.reclaim()
	    @everywhere GC.gc()
	end
        next!(prog)
	# Save the dictionaries to disk
	serialize(joinpath(folder, "efizz_coeffs_dict.jls"), coeffs_dict)
	serialize(joinpath(folder, "efizz_models_dict.jls"), models_dict)
	serialize(joinpath(folder, "efizz_inds_dict.jls"), inds_dict)
    end
end




p=plot(chns)
# Assuming you have 33 rows and 2 columns
for i in 1:33
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
