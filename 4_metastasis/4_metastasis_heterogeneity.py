# %% Import packages

import pandas as pd
import numpy as np
from pathlib import Path
import os
import plotly.express as px
import plotly.graph_objects as go

os.chdir(Path(__file__).parent)
output_dir = Path('.')
output_dir.mkdir(parents=True, exist_ok=True)

patients = ["PC06", "PC07", "PC08", "PC09", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16"]
condor_downstream_dir = Path("../0_condor_pipeline/condor_downstream/ete_trees_refined_subclonal_snvs")

# ----- read in sample anatomical name mapping file -----
sample_name_organ_map = pd.read_excel(
    '../Tapestri_batch2_samples_MASTER.xlsx', sheet_name="all_sample_clinical_info"
    ).set_index('sample')['HZ_official_site_name'].str.lower().str.replace(' ', '_').to_dict()

sample_name_general_site_map = pd.read_excel(
    '../Tapestri_batch2_samples_MASTER.xlsx', sheet_name="all_sample_clinical_info"
    ).set_index('sample')['Annotated Anatomic Site'].str.lower().str.replace(' ', '_').to_dict()
# series = series.str.lower().str.replace(' ', '_')

# %% Question 1: for each organ site, get: (1) how heterogeneous it is (2) how far it is from its corresponding primary tumor
def __get_heterogeneity_between_two_sites(site1_clone_compo, site2_clone_compo, clone_cn_profiles):
    """
    Calculate the heterogeneity between two sites, given their clone compositions and clone copy number profiles.
    site1_clone_compo: pandas series, the clone composition of site1
    site2_clone_compo: pandas series, the clone composition of site2
    clone_cn_profiles: pandas dataframe, the clone copy number profiles of all the clones in the patient
    """
    if not len(site1_clone_compo) == len(site2_clone_compo):
        raise ValueError("The length of the two clone compos are not equal!")
    # first calculate inter-clone distances
    # iterate through every pair of clones (rows) in clone_cn_profiles, get absolute difference
    site1_clone_compo.index = site1_clone_compo.index.astype(int)
    site2_clone_compo.index = site2_clone_compo.index.astype(int)
    clone_cn_profiles.index = clone_cn_profiles.index.astype(int)
    inter_cn_clone_dist = pd.DataFrame(columns=site1_clone_compo.index, index=site1_clone_compo.index)
    for clone_i in site1_clone_compo.index:
        for clone_j in site1_clone_compo.index:
            inter_cn_clone_dist.loc[clone_i, clone_j] = np.linalg.norm(clone_cn_profiles.loc[clone_i] - clone_cn_profiles.loc[clone_j])

    # distance between two sites is P1' x inter_cn_clone_dist x P2
    # P1' is the clone composition of site1, P2 is the clone composition of site2
    return np.dot(np.dot(site1_clone_compo, inter_cn_clone_dist), site2_clone_compo)

# %% Proceed patient by patient

site_stat_df = pd.DataFrame(columns=['patient', 'site_general','heterogeneity_within_self', 'heterogeneity_from_primary'])

pts_wo_primary = set()
for patient_i in patients:
    # ----- load clone compo -----
    clone_compo = pd.read_csv(condor_downstream_dir / patient_i / f"{patient_i}_clone_compo.csv")
    clone_compo = clone_compo.set_index("sample_name")
    # clone_compo.index = clone_compo.index.map(sample_name_map)
    clone_compo['site_general'] = clone_compo.index.map(sample_name_general_site_map)
    # clone_compo['site_specific'] = clone_compo.index.map(sample_name_organ_map)
    # collapse all the primary tumor samples into one
    clone_compo = clone_compo.groupby(['site_general']).sum()
    # normalize each clone (column) to sum to 1 
    clone_compo = clone_compo.div(clone_compo.sum(axis=1), axis=0)

    # ----- load condor clone cn_profiles -----
    condor_clone_cn_profiles = pd.read_csv(condor_downstream_dir / patient_i / f"{patient_i}_final_clone_cn_profiles.csv")
    condor_clone_cn_profiles.set_index(condor_clone_cn_profiles.columns[0], inplace=True)
    condor_clone_cn_profiles.index.name = 'clone_id'
    # # identify diploid clone:
    # diploid_clone_idx = condor_clone_cn_profiles[(condor_clone_cn_profiles== 2).all(axis=1)].index[0]

    # # remove the diploid clone column from clone_compo
    # clone_compo = clone_compo.drop(str(diploid_clone_idx), axis=1)

    patient_site_hetero_df = pd.DataFrame(index = clone_compo.index, columns = clone_compo.index)

    for site_i in clone_compo.index:
        for site_j in clone_compo.index:
            patient_site_hetero_df.loc[site_i, site_j] = __get_heterogeneity_between_two_sites(clone_compo.loc[site_i], clone_compo.loc[site_j], condor_clone_cn_profiles)

    # add to the site_stat_df
    for site_i in clone_compo.index:
        
        try: 
            dist_from_primary = patient_site_hetero_df.loc[site_i, 'primary']
        except KeyError:
            print(f"For patient {patient_i}, there is no primary tumor sample!")
            dist_from_primary = np.nan
            pts_wo_primary.add(patient_i)
        new_row = {'patient': patient_i, 'site_general': site_i, 'heterogeneity_within_self': patient_site_hetero_df.loc[site_i, site_i], 'heterogeneity_from_primary': dist_from_primary}
        # since Pandas 2.0, pandas.DataFrame.append() has been removed, use pandas.concat() instead
        site_stat_df = pd.concat([site_stat_df, pd.DataFrame(new_row, index=[0])], ignore_index=True)

    # if patient_i == 'RA15_16':
    #     break

# %% Mark the outlier
site_stat_df['outlier'] = False
site_stat_df.loc[site_stat_df['patient'] == 'PC10', 'outlier'] = True

# %% ===== Plot stripplots of clonal heterogeneity ===== 
import matplotlib.pyplot as plt
import seaborn as sns
from tea.plots import mpl_config_params

sites_of_interest = ['liver', 'peritoneum', 'lung', 'primary', 'diaphragm', 'lymph_node']

# plot stripplots (1x2)
mpl_config_params(font_size = 12)
fig, ax = plt.subplots(1, 2, figsize=(12, 8))
# A. heterogeneity within self
sns.boxplot(
    x='site_general', 
    y='heterogeneity_within_self', 
    # hue = 'outlier', palette = {False: 'royalblue', True: 'grey'},
    color="grey", boxprops=dict(alpha=.3), fliersize=0,
    data=site_stat_df[site_stat_df['site_general'].isin(sites_of_interest)], 
    ax=ax[0],
    order=sites_of_interest,
    # color = 'royalblue',
    # s=6, linewidth=1, edgecolor='black', alpha=.8,
    )
sns.stripplot(
    x='site_general', 
    y='heterogeneity_within_self', 
    # hue = 'outlier', palette = {False: 'royalblue', True: 'grey'},
    color = 'royalblue',
    data=site_stat_df[site_stat_df['site_general'].isin(sites_of_interest)], 
    ax=ax[0],
    order=sites_of_interest,
    # color = 'royalblue',
    s=6, linewidth=1, edgecolor='black', alpha=.8,
    )

# subtitle
ax[0].set_title("Heterogeneity within self", fontsize=16)
ax[0].set_xlabel('')
ax[0].set_ylabel('')
ax[0].set_ylim(0, 60)
# include the number of observations of each category in the xaxis labels
new_x_ticklabels = []
for tick in ax[0].get_xticklabels():
    tick.set_rotation(45)
    # tick.set_fontsize(12)
    # get the number of observations
    num_obs = len(site_stat_df[site_stat_df['site_general'] == tick.get_text()])
    tick.set_text(tick.get_text() + "\n" + f"n={num_obs}")
    new_x_ticklabels.append(tick.get_text())
ax[0].set_xticklabels(new_x_ticklabels)


# B. heterogeneity from primary
sns.boxplot(
    x='site_general', 
    y='heterogeneity_from_primary', 
    # hue = 'outlier', palette = {False: 'royalblue', True: 'grey'},
    color="grey", boxprops=dict(alpha=.3), fliersize=0,
    data=site_stat_df[site_stat_df['site_general'].isin(sites_of_interest)], 
    ax=ax[1],
    order=sites_of_interest,
    # color = 'royalblue',
    # s=6, linewidth=1, edgecolor='black', alpha=.8,
    )
sns.stripplot(
    x='site_general', 
    y='heterogeneity_from_primary', 
    hue = 'outlier', palette = {False: '#ffcc00', True: 'red'},
    data=site_stat_df[site_stat_df['site_general'].isin(sites_of_interest)], 
    ax=ax[1],
    order=sites_of_interest,
    color = 'orange',
    s=6, linewidth=1, edgecolor='black', alpha=.8,
    )
ax[1].set_title("Heterogeneity from primary", fontsize=16)
ax[1].set_xlabel('')
ax[1].set_ylabel('')
ax[1].set_ylim(0, 60)
# dsiable legend for hue
# ax[1].get_legend().remove()

# include the number of observations of each category in the xaxis labels
new_x_ticklabels = []
for tick in ax[1].get_xticklabels():
    tick.set_rotation(45)
    # tick.set_fontsize(12)
    # get the number of observations
    num_obs = len(site_stat_df[
        (site_stat_df["site_general"] == tick.get_text()) & 
        (site_stat_df["heterogeneity_from_primary"].notna())
	])
    tick.set_text(tick.get_text() + "\n" + f"n={num_obs}")
    new_x_ticklabels.append(tick.get_text())
ax[1].set_xticklabels(new_x_ticklabels)

sns.move_legend(ax[1], bbox_to_anchor=(1, 1.02), loc='upper left', title='PC10')

# fig.suptitle("Heterogeneity between sites", fontsize=20)
fig.tight_layout()
fig.subplots_adjust(top=0.85)
fig.savefig(output_dir / "organ_tropism_stripplots.pdf")
# %% ===== Draw pie charts =====

patient_names = ["PC09", "PC10", "PC11"]
for patient_name in patient_names:
    wd =  Path("sample_pie_charts")
    wd.mkdir(exist_ok=True, parents=True)
    clone_compo = condor_downstream_dir / patient_name / f"{patient_name}_clone_compo.csv"
    clone_compo_df = pd.read_csv(clone_compo, index_col=0)
    unique_cluster_ids_sorted = np.sort(clone_compo_df.columns.astype(int))

    cn_clone_palette = dict(zip(unique_cluster_ids_sorted, np.array(px.colors.qualitative.Pastel)[unique_cluster_ids_sorted]))

    # # %% Make a pie chart for each sample

    for sample_i in clone_compo_df.index.unique():
        sample_name_mapped = sample_name_organ_map[sample_i]
        # filter columns that have 0
        filtered_clone_compo_df = clone_compo_df.loc[sample_i].loc[clone_compo_df.loc[sample_i] > 0]
        fig = go.Figure(data=[go.Pie(
            labels=filtered_clone_compo_df.index, 
            values=filtered_clone_compo_df.values, 
            hole=.3
        )])
        colors = [cn_clone_palette[int(i)] for i in filtered_clone_compo_df.index]
        fig.update_traces(
            # hoverinfo='label+percent', 
            # textinfo='value', 
            # textinfo = 'label+percent',
            textfont_size=20,
            marker=dict(colors=colors, 
            line=dict(color='#000000', width=2))
        )
        fig.update_layout(
            font=dict(
                family="Arial",
                size=12,
            ),
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            showlegend=False,
        )
        # add title
        fig.add_annotation(
            text=f"{sample_name_mapped}",
            x=0.5,
            y=-0.2,
            showarrow=False,
            font=dict(
                size=40,
                color="black"
            )
        )
        fig.write_image(str(wd /f"{patient_name}_{sample_name_mapped}_clone_compo_pie.png"), scale=2)
# %%
