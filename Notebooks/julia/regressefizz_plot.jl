
plotfolder=joinpath(splitpath(folder)[1:end-1]..., "julia_regressefizz")
if !isdir(plotfolder)
	mkdir(plotfolder)
end

# PLOT: HEATMAPS: Ru components
#
function generate_heatmap(df::DataFrame, metric::Symbol;
    prop=:commsub,
    type=nothing, diff_type=nothing, kwargs...)
    
    # If diff_type is provided, ensure type is also provided
    if diff_type !== nothing && type === nothing
        throw(ArgumentError("If diff_type is provided, type must also be specified"))
    end
    
    # Create matrix for the heatmap
    matrix_data = Matrix{Float64}(undef, length(unique(df.Animal)), length(unique(df[!,prop])))

    df = sort(df, [:Animal, prop])
    
    for (i, animal) in enumerate(unique(df.Animal))
        for (j, behavior) in enumerate(unique(df[!,prop]))
            
            # If diff_type is not provided, use the original type behavior
            if diff_type === nothing
                subset = type === nothing ? df[(df.Animal .== animal) .& (df[!,prop] .== behavior), :] : 
                                             filter(row -> row.Type == type, df[(df.Animal .== animal) .& (df[!,prop] .== behavior), :])
                matrix_data[i, j] = mean(subset[:, metric])
                
            # If diff_type is provided, compute the difference
            else
                subset_type = filter(row -> row.Type == type, df[(df.Animal .== animal) .& (df[!,prop] .== behavior), :])
                subset_diff_type = filter(row -> row.Type == diff_type, df[(df.Animal .== animal) .& (df[!,prop] .== behavior), :])
                
                matrix_data[i, j] = mean(subset_type[:, metric]) - mean(subset_diff_type[:, metric])
            end
        end
    end
    
    title_str = diff_type === nothing ? "Heatmap of Mean $metric" : "Heatmap of Difference in $metric between $type and $diff_type"
    
    defaults = (;xlabel="Animal", ylabel="Commsub", title=title_str, color=:viridis, aspect_ratio=:auto)
    kw = merge(defaults, kwargs)
    
    # Generate heatmap
    heatmap(unique(df.Animal), unique(df[!,prop]), matrix_data'; kw...)
end

generate_heatmap(stat_df, :Accuracy, color=:balance; prop=:commsub)
savefig(joinpath(plotfolder, "heatmap_Accuracy_commsub.png"))
savefig(joinpath(plotfolder, "heatmap_Accuracy_commsub.pdf"))

generate_heatmap(stat_df, :R2, color=:balance, clim=(-1,1))
savefig(joinpath(plotfolder, "heatmap_R2.png"))
savefig(joinpath(plotfolder, "heatmap_R2.pdf"))

sort!(stat_df, [:Type,:Animal,:commsub])
plot_dict = OrderedDict()
type = first(unique(stat_df.Type))
for type in unique(stat_df.Type)
    plot_dict[type, :Accuracy] =
        generate_heatmap(stat_df, :Accuracy, color=:balance; type, title=string("Accuracy: $type"))
    plot_dict[type, :R2] =
        generate_heatmap(stat_df, :R2, color=:balance; type, title=string("R2: $type"))
end

# ----------------

summ = combine(groupby(stat_df, [:commsub, :Type]),
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


function generate_barplot(df::DataFrame; prop=:commsub, colors=nothing)
    # Create a grouped bar plot for each metric
    p1 = groupedbar(df[!,prop], df[:, :Accuracy_mean], group=df.Type, 
                    yerr=(df.Accuracy_mean .- df.Accuracy_lower, df.Accuracy_upper .- df.Accuracy_mean), 
                    title="Accuracy", legend=:topright, color_palette=colors)
hline!([0], color=:black, linestyle=:dash, label="")
    p2 = groupedbar(df[!,prop], df[:, :R2_mean], group=df.Type, 
                    yerr=(df.R2_mean .- df.R2_lower, df.R2_upper.-df.R2_mean), 
		title="R2", legend=:topright, ylim=(-0.5,0.75), color_palette=colors)
hline!([0], color=:black, linestyle=:dash, label="")
    # Combine the subplots
plot(p1, p2, layout=(2, 1), size=0.8.*(900,800), bbox=(0.2, 0.8), legend=:outertopright)
end

Base.adjoint(x::String) = x;
colors=get.([ColorSchemes.coolwarm], LinRange(0,1,length(unique(summ.Type)))|> collect)
p=generate_barplot(subset(summ, :commsub=>x->startswith.(x, "R") .&& vec(any(endswith.(x, ["1","2","3"]'), dims=2))  ), colors=colors)
savefig(joinpath(plotfolder, "barplot_R_commsub_Accuracy_R2_mean_comps=123.png"))
savefig(joinpath(plotfolder, "barplot_R_commsub_Accuracy_R2_mean_comps=123.pdf"))

p=generate_barplot(subset(summ, :commsub=>x->startswith.(x, "Ru") .&& vec(any(endswith.(x, ["1","2","3"]'), dims=2))  ), colors=colors)
savefig(joinpath(plotfolder, "barplot_Ru_commsub_Accuracy_R2_mean_comps=123.png"))
savefig(joinpath(plotfolder, "barplot_Ru_commsub_Accuracy_R2_mean_comps=123.pdf"))

# ----------------


