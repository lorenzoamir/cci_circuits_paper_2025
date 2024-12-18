import os
import re
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
import matplotlib.patheffects as path_effects  # Import PathEffects for text outlines

ncolor = 'C0'
tcolor = 'C1'

df = pd.read_csv('/home/lnemati/alya_survival/aggregate/better_pairs.csv', index_col=0)

def survival(interactions, df):
    import pandas as pd
    from lifelines import KaplanMeierFitter
    import re

    # Extract genes from the interaction string
    genes = list(set(re.split(r'[+&_]', interactions)))

    # Check if there are multiple conditions
    multiple_conditions = len(df['condition'].unique()) > 1

    # Filter and process the DataFrame
    cols = ['OS.time', 'OS', 'condition'] + genes
    df = df[cols]

    # Convert columns to numeric, except for condition
    df.loc[:, [col for col in df.columns if col != 'condition']] = (
        df.loc[:, [col for col in df.columns if col != 'condition']].apply(pd.to_numeric, errors='coerce')
    )
    df = df.dropna()

    # Convert OS time to years
    df['OS.time'] = df['OS.time'] / 365

    high_expression_group = pd.DataFrame()
    low_expression_group = pd.DataFrame()

    if multiple_conditions:
        # If there are multiple conditions, split the patients based on each condition
        for condition in df['condition'].unique():
            condition_df = df[df['condition'] == condition]

            # Calculate medians for each gene within the condition
            medians = condition_df[genes].median()

            # Identify patients above/below the median for all genes
            high_expression = condition_df[(condition_df[genes] > medians).T.all()]
            low_expression = condition_df[(condition_df[genes] <= medians).T.all()]

            # Append to the respective groups
            high_expression_group = pd.concat([high_expression_group, high_expression])
            low_expression_group = pd.concat([low_expression_group, low_expression])
    else:
        # If there is only one condition, split based on the entire dataset
        medians = df[genes].median()
        high_expression_group = df[(df[genes] > medians).T.all()]
        low_expression_group = df[(df[genes] <= medians).T.all()]

    # Count number of patients
    n_patients_low = low_expression_group.shape[0]
    n_patients_high = high_expression_group.shape[0]

    df['group'] = None
    df.loc[low_expression_group.index, 'group'] = 0
    df.loc[high_expression_group.index, 'group'] = 1
    df = df.dropna()

    # Fit the Kaplan-Meier estimators
    kmf_low = KaplanMeierFitter()
    kmf_high = KaplanMeierFitter()

    kmf_low.fit(
        low_expression_group['OS.time'],
        event_observed=low_expression_group['OS'],
        label="Low Expression"
    )
    kmf_high.fit(
        high_expression_group['OS.time'],
        event_observed=high_expression_group['OS'],
        label="High Expression"
    )

    return {
        'name': interactions,
        'kmf_low': kmf_low,
        'kmf_high': kmf_high,
        'low_expression': low_expression_group.index,
        'high_expression': high_expression_group.index,
    }

def get_genes(string):
    genes = list(re.split(r'[+&_]', string))

    # Remove duplicates while preserving order
    seen = set()
    ordered_genes = []
    for gene in genes:
        if gene not in seen:
            seen.add(gene)
            ordered_genes.append(gene)

    return ordered_genes

def plot_all_3(interaction, survival_df, savepath=None):
    if survival_df.tissue.nunique() == 1:
        tissue = survival_df.tissue[0].title().replace(' ', '_')
    else:
        tissue = 'Pan_Cancer'
    
    fig, axes = plt.subplots(3, 2, figsize=(12, 12))
    plt.subplots_adjust(hspace=0.5)
    interaction_results = survival(interaction, survival_df)
    
    int1, int2 = interaction.split('_')
    interaction_results1 = survival(int1, survival_df)
    interaction_results2 = survival(int2, survival_df)

    fs = 13
    
    # Plot the KM curves for both groups on the same plot
    for ax, result in zip(axes, (interaction_results, interaction_results1, interaction_results2)):
        name = result['name']
        genes = get_genes(name)
        
        pval = rawdf.set_index(['pair', 'tissue']).loc[(name, tissue), 'logrank_pval']
        n_low = rawdf.set_index(['pair', 'tissue']).loc[(name, tissue), 'n_patients_low']
        n_high = rawdf.set_index(['pair', 'tissue']).loc[(name, tissue), 'n_patients_high']
        #pval = result['p_value']
                
        ci = False
        censors=False
        result['kmf_low'].plot_survival_function(
            ax=ax[1],
            color=ncolor,
            ci_show=ci,
            show_censors=censors,
            label='Low Expression: {}'.format(n_low)
        )
        result['kmf_high'].plot_survival_function(
            ax=ax[1],
            color=tcolor,
            ci_show=ci,
            show_censors=censors,
            label='High Expression: {}'.format(n_high)

        )
        add_at_risk_counts(result['kmf_low'], result['kmf_high'], ax=ax[1], ypos=-0.35, rows_to_show=['At risk'])
        
        ax[0].set_title(name, fontsize=fs)
        ax[1].set_title(name, fontsize=fs)
        
        ax[1].set_xlabel('Time (Years)', fontsize=fs)
        ax[1].set_ylabel('Survival Prob.', fontsize=fs)
        #ax[1].set_xlim(0, 15)
        ax[1].set_ylim(0., 1.01)
        ax[1].legend(fontsize=fs)
        
        ypos = ax[1].get_ylim()[0] + 0.5*(ax[1].get_ylim()[1] - ax[1].get_ylim()[0])
        xpos = ax[1].get_xlim()[0] + 0.7*(ax[1].get_xlim()[1] - ax[1].get_xlim()[0])

        text = ax[1].text(xpos, ypos, '$p={:.1e}$'.format(pval), fontsize=fs)
        text.set_path_effects([
            path_effects.Stroke(linewidth=8, foreground=(1,1,1,0.5)),  # White outline
            path_effects.Normal()  # Normal text
        ])
        
        survival_df['Group'] = None
        survival_df.loc[survival_df.index.isin(result['high_expression']), 'Group'] = 'High Expression'
        survival_df.loc[survival_df.index.isin(result['low_expression']), 'Group'] = 'Low Expression'
        violindata = survival_df.melt(id_vars=['Group'], value_vars=genes)

        # Reorder violindata based on the genes list
        violindata['variable'] = pd.Categorical(violindata['variable'], categories=genes, ordered=True)

        # Create violin plot with explicit order
        sns.violinplot(
            data=violindata,
            x='variable',
            y='value',
            hue='Group',
            ax=ax[0],
            palette={'High Expression': tcolor, 'Low Expression': ncolor},
            order=genes,
            inner=None,
            cut=0
        )

        # Define the offset for the median line
        x_offset = 0.4
        lw = 2  # linewidth
        color = 'k'  # line color
        ls = '-'  # line style

        # Compute medians, preserving the order in genes
        median_values = (
            violindata.groupby('variable')['value'].median()
            .reindex(genes)  # Reorder by genes list
            .reset_index()
        )

        # Loop to add median lines with offset and collect coordinates
        #for i, (gene, median) in enumerate(zip(median_values['variable'], median_values['value'])):
        #    # Collect coordinates for the horizontal line (hlines)
        #        ax[0].plot([i - x_offset, i + x_offset], [median, median], color=color, linestyle=ls, linewidth=lw)

        #ax[0].set_xlim(-x_offset, len(genes)-1+x_offset)
        #ax[0].legend([])
        ax[0].set_xlabel('Gene', fontsize=fs)
        ax[0].set_ylabel('log2(0.001 + TPM)', fontsize=fs)
        ax[0].get_legend().remove()
        
        ax[0].tick_params(labelsize=fs)
        ax[1].tick_params(labelsize=fs)
        
        ax[0].spines[['bottom', 'left']].set_color('black')
        ax[0].spines[['right', 'top']].set_visible(False)
        ax[1].spines[['bottom', 'left']].set_color('black')
        ax[1].spines[['right', 'top']].set_visible(False)

    if savepath:
        plt.savefig(savepath, bbox_inches='tight')
    #plt.savefig('/home/lnemati/pathway_crosstalk/results/figures/survival/survival_example.pdf', bbox_inches='tight')
    
    plt.show()

rawdf = pd.read_csv('/home/lnemati/alya_survival/aggregate/all_unfiltered.csv')

# For each tissue, plot all
for tissue in df['tissue'].unique():
    tissue = tissue.replace(' ', '_')
    print(tissue)
    savedir = os.path.join('/home/lnemati/alya_survival/plots', tissue)
    os.makedirs(savedir, exist_ok=True)
    survival_df = pd.read_csv(os.path.join('/home/lnemati/pathway_crosstalk/data/survival_data/', tissue.lower() + '.csv'), index_col=0)
    for pair in df[df['tissue'] == tissue].index:
        print(pair)
        plot_all_3(pair, survival_df, savepath=os.path.join(savedir, pair + '.pdf'))

print('Done: plot.py')

