# assumes variables loaded from regressbeh.jl

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


# Coefficients by animal and behavior

using PyCall
ENV["MPLBACKEND"] = "Agg"
@pyimport seaborn as sns
@pyimport pandas as pd
@pyimport matplotlib.pyplot as plt
@pyimport mplcyberpunk
plt.style.use("cyberpunk")
# Convert the DataFrame to a dictionary
Cdf = Dict(names(coefs_df) .=> eachcol(coefs_df))
Cdf = pd.DataFrame(Cdf)
Cdf.to_csv(joinpath(folder, "coefs_df.csv"))
py"""
$(Cdf).query('Coefficient != "Intercept"', inplace=True)
"""
Cdf

# sns.set_style("cyberpunk")
sns.set_context("paper", font_scale=1.5, rc=Dict("lines.linewidth"=>2.5))
# sns.set_palette("colorblind")


# PLOT: COEFFICIENTS BY ANIMAL AND BEHAVIOR (BOXPLOTS), y=VALUE
g=sns.catplot(x="Coefficient", y="Value", hue="Animal", col="Behavior", data=Cdf, kind="box", height=4, aspect=2, legend_out=false, col_wrap=3)
g.set_xticklabels(rotation=45)
# Set ylims
for ax in g.axes
	ax.set_ylim(-1, 1)
	ax.axhline(0, color="black", linestyle="--")
	# Place a vertical line separating R and Ru labels
	ax.axvline(11.5, color="black", linestyle="--")
	# Remove bar edge color
	for bar in ax.patches
		bar.set_edgecolor("none")
	end
end

# PLOT: COEFFICIENTS BY ANIMAL AND BEHAVIOR (BOXPLOTS), y=ABSOLUTE VALUE
g=sns.catplot(x="Coefficient", y="ValueAbs", hue="Animal", col="Behavior", data=Cdf, kind="box", height=4, aspect=2, legend_out=false, col_wrap=3)
g.set_xticklabels(rotation=45)
# Set ylims
for ax in g.axes
	ax.set_ylim(0, 1)
	ax.axhline(0, color="black", linestyle="--")
	# Place a vertical line separating R and Ru labels
	ax.axvline(11.5, color="black", linestyle="--")
	# Remove bar edge color
	for bar in ax.patches
		bar.set_edgecolor("none")
	end
end
g.savefig(joinpath(plotfolder, "coefficients_boxplots_by_animal_behavior_y=abs.png"))
g.savefig(joinpath(plotfolder, "coefficients_boxplots_by_animal_behavior_y=abs.svg"))

g=sns.catplot(x="Coefficient", y="ValueAbs", col="Behavior", data=Cdf, kind="box", height=4, aspect=2, legend_out=false, col_wrap=3)
g.set_xticklabels(rotation=45)
# Set ylims
for ax in g.axes
	ax.set_ylim(0, 1)
	ax.axhline(0, color="black", linestyle="--")
	# Place a vertical line separating R and Ru labels
	ax.axvline(11.5, color="black", linestyle="--")
	# Remove bar edge color
	for bar in ax.patches
		bar.set_edgecolor("none")
	end
end

# import library reload 
py"""
import importlib
import matplotlib.pyplot as plt
importlib.reload(plt)
"""

g=sns.relplot(x="Coefficient", y="ValueAbs", col="Behavior", data=Cdf, kind="line", height=4, aspect=2, col_wrap=3, hue="orthostate", errorbar="ci", sharey=false)
g.set_xticklabels(rotation=45)
# g.map(sns.stripplot, x="Coefficient", y="ValueAbs", hue="Animal", data=Cdf, dodge=true, alpha=0.5, jitter=0.2)
# Set ylims
g.map(x -> mplcyberpunk.add_glow_effects(gradient_fill=true))
for ax in g.axes
	# ax.set_ylim(0, 1)
	ax.axhline(0, color="black", linestyle="--")
	# Place a vertical line separating R and Ru labels
	# ax.axvline(11.5, color="black", linestyle="--")
	# Remove bar edge color
	for bar in ax.patches
		bar.set_edgecolor("none")
	end
end
plt.show()

