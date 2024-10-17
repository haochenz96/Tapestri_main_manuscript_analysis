# %% [markdown]
# # Draw pie charts

# %%
cn_clone_assignment_df_updated_for_pie = cn_clone_assignment_df_updated.set_index('sample')

# %%
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path


# load clone assignment df
condor_downstream_dir = Path("/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/0_condor_pipeline/condor_downstream")

# load MASTER df
sample_sheet = Path("../Tapestri_batch2_samples_MASTER_INTERNAL.xlsx")
sample_sheet_df = pd.read_excel(sample_sheet, sheet_name="all_sample_clinical_info")
sample_name_organ_map = dict(zip(sample_sheet_df['sample'], sample_sheet_df['HZ_official_site_name']))

# %% Load clone assignment df

patient_names = ["RA16_29", "RA17_13", "RA17_22", "RA21_17"]
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
        fig.write_image(str(wd /f"{patient_name}_{sample_name_mapped}_clone_compo_pie.png"), scale=2)# %%

# %%
