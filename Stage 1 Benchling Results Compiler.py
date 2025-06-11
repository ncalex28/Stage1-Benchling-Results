#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px


#### Functions

def merge_and_clean(df1, df2, keys):
    """Merge two DataFrames on given keys, drop duplicate '_df2' columns, and clean '_df1' suffixes."""
    df_merged = pd.merge(df1, df2, on=keys, how="outer", suffixes=("_df1", "_df2"))
    df_merged = df_merged.drop(columns=[col for col in df_merged.columns if "_df2" in col])
    df_merged = df_merged.rename(columns=lambda x: x.replace("_df1", ""))
    return df_merged


### Input file paths

# Use latest export from Benchling
collection_path= r"C:\Users\nalexander\Desktop\Python\2025\Stage1 Export\20250610\2025-06-10 03_58_00 PM PD Stage1 Dashboard - Viable Cell Density (cells-ml).csv"
metaflex_path= r"C:\Users\nalexander\Desktop\Python\2025\Stage1 Export\20250610\2025-06-10 03_58_23 PM PD Stage1 Dashboard - Metaflex_ Lactate -time.csv"
flow_path= r"C:\Users\nalexander\Desktop\Python\2025\Stage1 Export\20250610\2025-06-10 03_58_13 PM PD Stage1 Dashboard - CD34+-KDR+(%).csv"
# agg_path= r"C:\Users\nalexander\Desktop\Python\2025\Stage1 Export\20250327\2025-03-27 09_39_17 AM PD Stage1 Dashboard - Mean Aggregate Size (um).csv"

# Keep this path constant
user_path= r"C:\Users\nalexander\Desktop\Python\2025\Stage1 Export\2025-03-27 10_51_24 AM PD Stage1 Dashboard - Experiment Creators.csv"


## Define .csv save path
save_path = r"C:\Users\nalexander\Desktop\Python\Benchling\Stage1 Exports\20250610_Stage1 Export_"


## Bulk code

## Import dfs
collection_df= pd.read_csv(collection_path)
metaflex_df= pd.read_csv(metaflex_path)
flow_df= pd.read_csv(flow_path)
# agg_df= pd.read_csv(agg_path)
user_df= pd.read_csv(user_path)

# Clean columns
mcols= ['EffectorID', 'Alias', 'Parent IPSC', 'Experiment ID', 'Day of Differentiation', 'Vessel Type', 'pH', 'Glu(g/L)', 'Lac(g/L)']
metaflex_df= metaflex_df.loc[:, mcols]

# Define merge keys
merge_keys = ["Experiment ID", "EffectorID", "Day of Differentiation"]

# Perform sequential merges for data tables
df_merge = collection_df
for df in [metaflex_df, flow_df, 
#            agg_df
          ]:
    df_merge = merge_and_clean(df_merge, df, merge_keys)

# Perform merge for experiment creator
df_merge = df_merge.merge(user_df, on= 'Experiment ID')


### Filter for only the Stage1 Experiments
##### These will presumably contain the text PPD-HPC-24XXXX

df_merge = df_merge.replace({'PPD-HSPC-25008': 'PPD-HPC-25008'}, regex=True)
df_merge = df_merge[df_merge["Experiment ID"].str.contains(r"PPD-HPC-24|PPD-HPC-25", na=False)]


## Melt
id_cols= ['EffectorID', 'Alias', 'Parent IPSC', 'Experiment ID', 'Experiment Created at', 'Experiment Creator', 'Vessel Type', 'Day of Differentiation', 'Culture Volume', 'Cytometer']
val_cols= ['Viability (%)', 'Viable Cell Density (cells/ml)', 'Total Cell Yield', 'Fold Expansion',
       'pH', 'Glu(g/L)', 'Lac(g/L)', 'Cytometer', 'CD34+(%)',
       'CD34+/CD43+(%)', 'CD34+/KDR+(%)', 'CD34+/CD31+(%)', 'CD34+/CXCR4+(%)',
       'CD31+/KDR+(%)', 'CD34+/CD62L+(%)', 'CD34+/Notch1+(%)',
       'CD34+/CD36+(%)', 'CD34+/CD41+(%)', 'CD36+(%)', 'CD41+(%)']

df_melt= pd.melt(df_merge, id_vars= id_cols, value_vars= val_cols)
df_melt= df_melt[df_melt['value'].notna()] #drop missing data
df_melt= df_melt[df_melt['variable'].notna()] #drop missing variables even though missing variables shouldn't be possible


## Pivot for correlations

# Ensure 'Day of Differentiation' is numeric for proper sorting (if necessary)
df_merge['Day of Differentiation'] = pd.to_numeric(df_merge['Day of Differentiation'], errors='coerce')

# Define metadata columns (all columns that are NOT measurements)
metadata_cols = ['EffectorID', 'Parent IPSC', 'Alias', 'Experiment ID', 'Vessel Type', 'Culture Volume', 'Experiment Creator']

# Define measurement columns (all columns that should be pivoted)
value_vars = [col for col in df_merge.columns if col not in metadata_cols + ['Day of Differentiation']]

# Reshape using melt (long format)
df_melted = df_merge.melt(id_vars=metadata_cols + ['Day of Differentiation'], 
                           value_vars=value_vars, 
                           var_name='Variable', 
                           value_name='Value')

# Filter only numeric columns for pivoting
df_melted['Value'] = pd.to_numeric(df_melted['Value'], errors='coerce')  # Convert values to numeric, forcing errors to NaN

# Pivot to wide format, filtering numeric columns
df_pivot = df_melted.pivot_table(index=metadata_cols, 
                                 columns=['Day of Differentiation', 'Variable'], 
                                 values='Value', 
                                 aggfunc='mean')  # Ensure aggregation for duplicates

# Flatten MultiIndex columns into a single level (with a proper naming convention)
df_pivot.columns = ['D{}_{}'.format(int(day), var) for day, var in df_pivot.columns]
df_pivot = df_pivot.reset_index()

# Show transformed DataFrame
df_pivot.head()


### Export dataframes to .csv
df_merge.to_csv(save_path+"_merged dataframe.csv")
df_melt.to_csv(save_path+"_melted dataframe.csv")
df_pivot.to_csv(save_path+"_pivot dataframe.csv")



#### Flexible filtering function for Streamlit heatmaps
## Flexible filter for df to flow heatmap

def plot_filtered_heatmap(df, id_column, filters, flow_columns, figsize=(12, 6), cmap='viridis'):
    """
    Takes merged df from Benchling Stage1 tables and plots heatmap of flow results
    Also returns the filtered df
    
    id_column= Effector ID, Alias, etc
    filters= Day of Differentiation, Experiment ID, vessel, etc
    cmap= viridis, mako, coolwarm, etc
    """
    
    #filtering
    filtered_df = df.copy()
    for key, value in filters.items():
        if callable(value):
            filtered_df = filtered_df[filtered_df[key].apply(value)]
        elif isinstance(value, list):
            filtered_df = filtered_df[filtered_df[key].isin(value)]
        else:
            filtered_df = filtered_df[filtered_df[key] == value]

    #set index and select flow columns
    filtered_df = filtered_df.set_index(id_column)
    heatmap_data = filtered_df[flow_columns]

    #plot heatmap
    plt.figure(figsize= figsize)
    sns.heatmap(heatmap_data, annot= True, cmap= cmap)
    plt.title(f"Filtered Heatmap ({len(heatmap_data)} samples)")
    plt.xlabel("Flow Markers")
    plt.ylabel(id_column)
    plt.tight_layout()
    plt.show()

    return filtered_df  #return df


# In[41]:


## Define filters
filters = {
    'Experiment ID': ["PPD-HPC-25003: Generation of CL308 Stage1 Demo Material from Branchburg Mock MCB", 
                      'PPD-HPC-25005: CNTY257 Demo 2 Run - Branchburg'],
    'Day of Differentiation': 14
}

## Define flow columns
flow_columns = ['CD34+(%)', 'CD34+/CD43+(%)',
       'CD34+/KDR+(%)', 'CD34+/CD31+(%)', 'CD34+/CXCR4+(%)', 'CD31+/KDR+(%)',
       'CD34+/CD62L+(%)', 'CD34+/Notch1+(%)', 'CD34+/CD36+(%)',
       'CD34+/CD41+(%)', 'CD36+(%)', 'CD41+(%)']


# Plot
plot_filtered_heatmap(
    df=df_merge,
    id_column='Alias',
    filters=filters,
    flow_columns=flow_columns,
    figsize=(12, 6),
    cmap='viridis'
)

print("Filters:")
for key, value in filters.items():
    print(f"{key}: {value}")





