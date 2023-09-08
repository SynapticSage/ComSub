
# ========
# ALERT:
# UV magnitude OVER TIME 󰥔
# ========

# TODO: CUT OUT DIFF GREATER THAN 1 minute

rhythms = df_clean.pattern_class_name.unique()
rhythm = rhythms[0]
for rhythm in rhythms:
    figfolder_lfp = os.path.join(figfolder, rhythm)
    if not os.path.exists(figfolder_lfp):
        os.makedirs(figfolder_lfp)
    df_matrix = prep_uv_melt(df_clean.query('pattern_class_name == @rhythm & pattern_cca2 == @overall_pattern'))
    df_matrix = prep_uv_magnitude_over_time(df_matrix)
    df_matrix.query('pattern_cca2 == @overall_pattern', inplace=True)
    df_matrix = df_matrix.reset_index()
    df_matrix['magnitude_u'] = np.sqrt(np.sum(df_matrix[[f"{i}.0_u" for i in range(1, 6)]], axis=1))
    df_matrix['magnitude_v'] = np.sqrt(np.sum(df_matrix[[f"{i}.0_v" for i in range(1, 6)]], axis=1))
    df_matrix['magnitude_r'] = np.sqrt(df_matrix['magnitude_u'] ** 2 + df_matrix['magnitude_v'] ** 2)
    df_matrix['on_r_commsub'] = (df_matrix['magnitude_u'] + df_matrix['magnitude_v']) / np.sqrt(2)
    df_matrix['off_r_commsub'] = (df_matrix['magnitude_u'] - df_matrix['magnitude_v']) / np.sqrt(2)
    df_matrix['on_r_commsub_div_total'] = df_matrix['on_r_commsub'] / df_matrix['magnitude_r']
    df_matrix['on_versus_off'] = np.log(df_matrix['on_r_commsub'] / df_matrix['off_r_commsub'])
    for i in range(1, 6):
        df_matrix[f"magnitude_r_{i}"] = df_matrix[f"{i}.0_u"] ** 2 + df_matrix[f"{i}.0_v"] ** 2
        df_matrix[f"on_r_commsub_{i}"] = (df_matrix[f"{i}.0_u"] + df_matrix[f"{i}.0_v"]) / np.sqrt(2)
        df_matrix[f"off_r_commsub_{i}"] = (df_matrix[f"{i}.0_u"] - df_matrix[f"{i}.0_v"]) / np.sqrt(2)
        df_matrix[f"on_r_commsub_div_total_{i}"] = df_matrix[f"on_r_commsub_{i}"] / df_matrix[f"magnitude_r_{i}"]
        df_matrix[f"on_versus_off_{i}"] = np.log(df_matrix[f"on_r_commsub_{i}"] / df_matrix[f"off_r_commsub_{i}"])

    df_matrix.dropna(subset=['on_r_commsub', 'off_r_commsub'], axis=0, inplace=True)
    df_matrix.query('epoch >= 2', inplace=True)

    # Set early and late epoch labels (first and last half)
    df_matrix['epoch'] = df_matrix['epoch'].astype(int)
    df_matrix['epoch_half'] = df_matrix['epoch'] % 4
    df_matrix['epoch_half'] = df_matrix['epoch_half'].replace({0: 'early', 2: 'late'})


    #9 Melt the data for easy plotting
    df_melted =
    df_matrix.melt(id_vars=['epoch','time','animal','genH','highlow','epoch_half'],
                   value_vars=['on_r_commsub','off_r_commsub',
                               'on_r_commsub_1', 'off_r_commsub_1',
                               'on_r_commsub_2', 'off_r_commsub_2',
                               'on_r_commsub_3', 'off_r_commsub_3',
                               'on_r_commsub_4', 'off_r_commsub_4',
                               'on_r_commsub_5', 'off_r_commsub_5'],
                   var_name='component', value_name='magnitude')
    min_animal_samples = df_melted.groupby(['animal']).size().min()
    df_melted = df_melted.groupby(['animal']).sample(min_animal_samples)
    # Plot using relplot
    df_melted.sort_values(['component','genH','highlow'], inplace=True)

    print("Plotting...")
    g = sns.relplot(data=df_melted, x='epoch', y='magnitude', kind='line',
                    hue='genH', col="highlow", height=3, aspect=0.5, row="component",
                    facet_kws={'sharey':False}, col_order=['low', 'mid', 'high'])
    for ax in g.axes.flat:
        # ax.axhline([0], color='black', linestyle='--')
        # replace | with \n
        ax.set_title(ax.get_title().replace('|', '\n'))
    # Equalize ylim on each row of axes
    for axs in g.axes:
        # Get the max ylim
        max_ylim = np.max([ax.get_ylim() for ax in axs])
        min_ylim = np.min([ax.get_ylim() for ax in axs])
        # Set all ylims to max_ylim
        for ax in axs:
            ax.set_ylim(min_ylim, max_ylim)
    g.set_axis_labels("Epoch", "Magnitude")
    g.tight_layout()
    plt.show()
    plt.savefig(os.path.join(figfolder_lfp,'match_to_genH_overtime_hue=genH_col=highlow_row=component.png'), dpi=300)
    plt.savefig(os.path.join(figfolder_lfp,'match_to_genH_overtime_hue=genH_col=highlow_row=component.pdf'))


    # Let's do a similar plot but with x=epoch_half
    g = sns.relplot(data=df_melted, x='epoch_half', y='magnitude', kind='line',
                    hue='genH', col="highlow", height=3, aspect=0.5, row="component",
                    facet_kws={'sharey':False}, col_order=['low', 'mid', 'high'],
                    markers='o')
    for ax in g.axes.flat:
        # ax.axhline([0], color='black', linestyle='--')
        # replace | with \n
        ax.set_title(ax.get_title().replace('|', '\n').replace('component = ',''))
        # make title size smaller
        ax.title.set_size(9)
        # make shading of ribbons darker
        for c in ax.collections:
            c.set_alpha(0.5)
    # Add more vspace between rows
    g.fig.subplots_adjust(hspace=1.2)
    # Equalize ylim on each row of axes
    for axs in g.axes:
        # Get the max ylim
        max_ylim = np.max([ax.get_ylim() for ax in axs])
        min_ylim = np.min([ax.get_ylim() for ax in axs])
        # Set all ylims to max_ylim
        for ax in axs:
            ax.set_ylim(min_ylim, max_ylim)
    g.set_axis_labels("Epoch", "Magnitude")
    g.tight_layout()
    plt.show()
    # Differentiate the save name from above
    plt.savefig(os.path.join(figfolder_lfp,'match_to_genH_epoch_half_overtime_hue=genH_col=highlow_row=component.png'), dpi=300)
    plt.savefig(os.path.join(figfolder_lfp,'match_to_genH_epoch_half_overtime_hue=genH_col=highlow_row=component.pdf'))

    # condensed version excluding individual components of the above

    df_melted = df_matrix.melt(id_vars=['epoch','time','animal','genH','highlow','epoch_half'], value_vars=['on_r_commsub','off_r_commsub'], var_name='component', value_name='magnitude')
    min_animal_samples = df_melted.groupby(['animal']).size().min()
    df_melted = df_melted.groupby(['animal']).sample(min_animal_samples)
    # Plot using relplot
    df_melted.sort_values(['component','genH','highlow'], inplace=True)

    print("Plotting...")
    g = sns.relplot(data=df_melted, x='epoch', y='magnitude', kind='line',
                    hue='genH', col="highlow", height=3, aspect=0.5, row="component",
                    facet_kws={'sharey':False}, col_order=['low', 'mid', 'high'])
    for ax in g.axes.flat:
        # replace | with \n
        ax.set_title(ax.get_title().replace('|', '\n'))
    # Equalize ylim on each row of axes
    for axs in g.axes:
        # Get the max ylim
        max_ylim = np.max([ax.get_ylim() for ax in axs])
        min_ylim = np.min([ax.get_ylim() for ax in axs])
        # Set all ylims to max_ylim
        for ax in axs:
            ax.set_ylim(min_ylim, max_ylim)
    g.set_axis_labels("Epoch", "Magnitude")
    g.tight_layout()
    plt.show()
    plt.savefig(os.path.join(figfolder_lfp,'condensed_match_to_genH_overtime_hue=genH_col=highlow_row=component.png'), dpi=300)
    plt.savefig(os.path.join(figfolder_lfp,'condensed_match_to_genH_overtime_hue=genH_col=highlow_row=component.pdf'))
    plt.close('all')

    # Let's do a similar plot but with x=epoch_half
    g = sns.relplot(data=df_melted, x='epoch_half', y='magnitude', kind='line',
                    hue='genH', col="highlow", height=3, aspect=0.5, row="component",
                    facet_kws={'sharey':False}, col_order=['low', 'mid', 'high'],
                    markers='o')
    for ax in g.axes.flat:
        # replace | with \n
        ax.set_title(ax.get_title().replace('|', '\n').replace('component = ',''))
        # make title size smaller
        ax.title.set_size(9)
        # make shading of ribbons darker
        for c in ax.collections:
            c.set_alpha(0.5)
    # Add more vspace between rows
    g.fig.subplots_adjust(hspace=1.2)
    # Equalize ylim on each row of axes
    for axs in g.axes:
        # Get the max ylim
        max_ylim = np.max([ax.get_ylim() for ax in axs])
        min_ylim = np.min([ax.get_ylim() for ax in axs])
        # Set all ylims to max_ylim
        for ax in axs:
            ax.set_ylim(min_ylim, max_ylim)
    g.set_axis_labels("Epoch", "Magnitude")
    g.tight_layout()
    plt.show()
    # Differentiate the save name from above
    plt.savefig(os.path.join(figfolder_lfp,'condensed_match_to_genH_epoch_half_overtime_hue=genH_col=highlow_row=component.png'), dpi=300)
    plt.savefig(os.path.join(figfolder_lfp,'condensed_match_to_genH_epoch_half_overtime_hue=genH_col=highlow_row=component.pdf'))
    plt.close('all')


    # Let's do a similar plot but with x=epoch_half
    g = sns.relplot(data=df_melted, x='epoch_half', y='magnitude', kind='line',
                    col='genH', hue="highlow", height=3, aspect=0.5, row="component",
                    facet_kws={'sharey':False}, hue_order=['low', 'mid', 'high'],
                    markers='o')
    for ax in g.axes.flat:
        # replace | with \n
        ax.set_title(ax.get_title().replace('|', '\n').replace('component = ',''))
        # make title size smaller
        ax.title.set_size(9)
        # make shading of ribbons darker
        for c in ax.collections:
            c.set_alpha(0.5)
    # Add more vspace between rows
    g.fig.subplots_adjust(hspace=1.2)
    # Equalize ylim on each row of axes
    for axs in g.axes:
        # Get the max ylim
        max_ylim = np.max([ax.get_ylim() for ax in axs])
        min_ylim = np.min([ax.get_ylim() for ax in axs])
        # Set all ylims to max_ylim
        for ax in axs:
            ax.set_ylim(min_ylim, max_ylim)
    g.set_axis_labels("Epoch", "Magnitude")
    g.tight_layout()
    plt.show()
    # Differentiate the save name from above
    plt.savefig(os.path.join(figfolder_lfp,'condensed_match_to_genH_epoch_half_overtime_hue=highlow_col=genH_row=component.png'), dpi=300)
    plt.savefig(os.path.join(figfolder_lfp,'condensed_match_to_genH_epoch_half_overtime_hue=highlow_col=genH_row=component.pdf'))
    plt.close('all')

    df_melted.sort_values(['highlow'], inplace=True, ascending=True, key=lambda x: x.map({'low': 0, 'mid': 1, 'high': 2}))
    g=sns.catplot(data=df_melted, x='highlow', y='magnitude', kind='box',
                  hue='genH', col='component', col_wrap=3, height=3, aspect=0.5,
                  hue_order=['coherence','power'])
    plt.savefig(os.path.join(figfolder_lfp,'catplot_condensed_match_to_genH_epoch_half_overtime_hue=genH_col=component.png'), dpi=300)
    plt.savefig(os.path.join(figfolder_lfp,'catplot_condensed_match_to_genH_epoch_half_overtime_hue=genH_col=component.pdf'))
    plt.close('all')




    # df_melted = df_matrix.melt(id_vars=['time','epoch','animal','genH'], value_vars=['on_r_comsub_smooth', 'off_r_comsub_smooth'], var_name='component', value_name='magnitude')
    # # Plot using relplot
    # print("Plotting...")
    # g = sns.relplot(data=df_melted, x='time', y='magnitude', kind='line',
    #                 hue='component', height=5, aspect=2, alpha=0.5, col='epoch',
    #                 row='animal', facet_kws=dict(sharex=False))
    # g.set_axis_labels("Event Time", "Magnitude")
    # g.tight_layout()
    # plt.show()
    # plt.savefig(os.path.join(figfolder_lfp,'match_to_genH_overtime_hue=component_row=animal_col=epoch.png'), dpi=300)
    # plt.savefig(os.path.join(figfolder_lfp,'match_to_genH_overtime_hue=component_row=animal_col=epoch.pdf'))
    #
    # g = sns.relplot(data=df_melted, x='time', y='magnitude', kind='line',
    #                 hue='genH', height=5, aspect=2, alpha=0.5, col='epoch',
    #                 row='animal', facet_kws=dict(sharex=False))
    # g.set_axis_labels("Event Time", "Magnitude")
    # g.tight_layout()
    # plt.show()
    # plt.savefig(os.path.join(figfolder_lfp,'match_to_genH_overtime_hue=component_row=animal_col=epoch.png'), dpi=300)
    # plt.savefig(os.path.join(figfolder_lfp,'match_to_genH_overtime_hue=component_row=animal_col=epoch.pdf'))
    # plt.close('all')
    #
    #
    # # Plot using relplot
    # print("Plotting...")
    # g = sns.relplot(data=df_matrix, x='time_bin', y='on_r_comsub', kind='line', hue='genH', height=5, aspect=2)
    # g.set_axis_labels("Time Bin", "Magnitude")
    # g.tight_layout()
    # plt.show()
    # # different save name from above and below with enough description to know what's being plotted
    # plt.savefig(os.path.join(figfolder_lfp,'match_to_genH_overtime_hue=genH.png'), dpi=300)
    # plt.savefig(os.path.join(figfolder_lfp,'match_to_genH_overtime_hue=genH.pdf'))
    #
    # # Plot using relplot
    # print("Plotting...")
    # g = sns.relplot(data=df_matrix, x='time_bin', y='on_r_comsub', kind='line', hue='genH_highlow', height=5, aspect=2)
    # g.set_axis_labels("Time Bin", "Magnitude")
    # g.tight_layout()
    # plt.show()
    # plt.savefig(os.path.join(figfolder_lfp,'match_to_genH_overtime_hue=genH_highlow.png'), dpi=300)
    # plt.savefig(os.path.join(figfolder_lfp,'match_to_genH_overtime_hue=genH_highlow.pdf'))
    #
    #
    #
    # # Plot using relplot
    # print("Plotting...")
    # g = sns.relplot(data=df_matrix, x='time_bin', y='off_r_comsub', kind='line', hue='genH', height=5, aspect=2)
    # g.set_axis_labels("Time Bin", "Magnitude")
    # g.tight_layout()
    # plt.show()
    #
    #
    # # Plot using relplot
    # print("Plotting...")
    # g = sns.relplot(data=df_matrix, x='time_bin', y='off_r_comsub', kind='line', hue='genH_highlow', height=5, aspect=2)
    # g.set_axis_labels("Time Bin", "Magnitude")
    # g.tight_layout()
    # plt.show()

    # Creat smoothed versions of our components over time_bin
    df_matrix.sort_values(['animal', 'genH', 'highlow', 'time'], inplace=True)
    df_matrix.loc[:,'on_r_comsub_smooth'] = df_matrix.groupby(['animal', 'genH', 'highlow'])['on_r_commsub'].transform(lambda x: x.rolling(200, 1).mean()).interpolate()
    df_matrix['off_r_comsub_smooth'] = df_matrix.groupby(['animal', 'genH', 'highlow'])['off_r_commsub'].transform(lambda x: x.rolling(200, 1).mean())
    # Create a smooth version for each of the subcomponents 1-5
    for i in range(1, 6):
        df_matrix[f"on_r_comsub_{i}_smooth"] = df_matrix.groupby(['animal', 'genH', 'highlow'])[f"on_r_commsub_{i}"].transform(lambda x: x.rolling(200, 1).mean())
        df_matrix[f"off_r_comsub_{i}_smooth"] = df_matrix.groupby(['animal', 'genH', 'highlow'])[f"off_r_commsub_{i}"].transform(lambda x: x.rolling(200, 1).mean())

    # Plot using relplot
    print("Plotting...")
    g = sns.relplot(data=df_matrix.query('genH == "coherence"'), x='time_bin', y='on_r_comsub_smooth', kind='line', hue='genH_highlow', height=5, aspect=2)
    g.set_axis_labels("Time Bin", "Magnitude")
    g.tight_layout()
    plt.show()
    # come up with descriptive name that does not follow the format above
    plt.savefig(os.path.join(figfolder_lfp,'time_bin_on_r_commsub_genH=coherence_overtime_hue=genH_highlow.png'), dpi=300)
    plt.savefig(os.path.join(figfolder_lfp,'time_bin_on_r_commsub_genH=coherence_overtime_hue=genH_highlow.pdf'))
           

    # now for off
    print("Plotting...")
    g = sns.relplot(data=df_matrix.query('genH == "coherence"'), x='time_bin', 
                    y='off_r_comsub_smooth', kind='line', hue='genH_highlow', height=5, aspect=2)
    g.set_axis_labels("Time Bin", "Magnitude")
    g.tight_layout()
    plt.show()
    # come up with descriptive name that does not follow the format above
    plt.savefig(os.path.join(figfolder_lfp,'time_bin_off_r_commsub_genH=coherence_overtime_hue=genH_highlow.png'), dpi=300)

    # THe above two paragraphs split by animal along the row
    # Plot using relplot
    print("Plotting...")
    g = sns.relplot(data=df_matrix.query('genH == "coherence"'), x='time_bin', y='on_r_comsub_smooth', kind='line', hue='genH_highlow', height=5, aspect=2, row='animal', facet_kws=dict(sharex=False))
    g.set_axis_labels("Time Bin", "Magnitude")
    g.tight_layout()
    plt.show()
    plt.savefig(os.path.join(figfolder_lfp,'time_bin_on_r_commsub_genH=coherence_overtime_hue=genH_highlow_row=animal.png'), dpi=300)
    plt.savefig(os.path.join(figfolder_lfp,'time_bin_on_r_commsub_genH=coherence_overtime_hue=genH_highlow_row=animal.pdf'))

    # now for off
    print("Plotting...")
    g = sns.relplot(data=df_matrix.query('genH == "coherence"'), x='time_bin', y='off_r_comsub_smooth', kind='line', hue='genH_highlow', height=5, aspect=2, row='animal', facet_kws=dict(sharex=False))
    g.set_axis_labels("Time Bin", "Magnitude")
    g.tight_layout()
    plt.show()
    plt.savefig(os.path.join(figfolder_lfp,'time_bin_off_r_commsub_genH=coherence_overtime_hue=genH_highlow_row=animal.png'), dpi=300)
    plt.savefig(os.path.join(figfolder_lfp,'time_bin_off_r_commsub_genH=coherence_overtime_hue=genH_highlow_row=animal.pdf'))

    plt.close('all')

    # Now plot on and off component i (animal summary, not each)
    for i in range(1, 6):
        print("Plotting...")
        g = sns.relplot(data=df_matrix.query('genH == "coherence"'), x='time_bin', y=f'on_r_comsub_{i}_smooth', kind='line', hue='genH_highlow', height=5, aspect=2)
        g.set_axis_labels("Time Bin", "Magnitude")
        g.tight_layout()
        plt.show()
        plt.savefig(os.path.join(figfolder_lfp,f'{i}_time_bin_on_r_commsub_genH=coherence_overtime_hue=genH_highlow.png'), dpi=300)
        plt.savefig(os.path.join(figfolder_lfp,f'{i}_time_bin_on_r_commsub_genH=coherence_overtime_hue=genH_highlow.pdf'))

        # now for off
        print("Plotting...")
        g = sns.relplot(data=df_matrix.query('genH == "coherence"'), x='time_bin', y=f'off_r_comsub_{i}_smooth', kind='line', hue='genH_highlow', height=5, aspect=2)
        g.set_axis_labels("Time Bin", "Magnitude")
        g.tight_layout()
        plt.show()
        plt.savefig(os.path.join(figfolder_lfp,f'{i}_time_bin_off_r_commsub_genH=coherence_overtime_hue=genH_highlow.png'), dpi=300)
        plt.savefig(os.path.join(figfolder_lfp,f'{i}_time_bin_off_r_commsub_genH=coherence_overtime_hue=genH_highlow.pdf'))
        plt.close('all')

    # And now with animals
    for i in range(1, 6):
        print("Plotting...")
        g = sns.relplot(data=df_matrix.query('genH == "coherence"'), x='time_bin', y=f'on_r_comsub_{i}_smooth', kind='line', hue='genH_highlow', height=5, aspect=2, row='animal', facet_kws=dict(sharex=False))
        g.set_axis_labels("Time Bin", "Magnitude")
        g.tight_layout()
        plt.show()
        plt.savefig(os.path.join(figfolder_lfp,f'{i}_time_bin_on_r_commsub_genH=coherence_overtime_hue=genH_highlow_row=animal.png'), dpi=300)
        plt.savefig(os.path.join(figfolder_lfp, f'{i}_time_bin_on_r_commsub_genH=coherence_overtime_hue=genH_highlow_row=animal.pdf'))

        # now for off
        print("Plotting...")
        g = sns.relplot(data=df_matrix.query('genH == "coherence"'), x='time_bin', y=f'off_r_comsub_{i}_smooth', kind='line', hue='genH_highlow', height=5, aspect=2, row='animal', facet_kws=dict(sharex=False))
        g.set_axis_labels("Time Bin", "Magnitude")
        g.tight_layout()
        plt.show()
        plt.savefig(os.path.join(figfolder_lfp, f'{i}_time_bin_off_r_commsub_genH=coherence_overtime_hue=genH_highlow_row=animal.png'), dpi=300)
        plt.savefig(os.path.join(figfolder_lfp, f'{i}_time_bin_off_r_commsub_genH=coherence_overtime_hue=genH_highlow_row=animal.pdf'))
        plt.close('all')
