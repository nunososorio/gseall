

import streamlit as st
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

# Define a function to perform GSEA and plot the results for a given database
def perform_gsea(database, genes, organism):
    # Run GSEA analysis
    result = gp.gsea(
        data = genes,
        gene_sets = database,
        organism = organism,
        permutation_num = 1000,
        method = 'signal_to_noise',
        processes = 1,
        seed = 1,
        outdir = None,
        no_plot = True
    )
    # Filter the significant terms
    significant_terms = result.results.loc[result.results['fdr'] < 0.05]
    # Plot a dot plot for the significant terms
    if not significant_terms.empty:
        fig, ax = plt.subplots()
        ax.scatter(
            significant_terms['es'],
            significant_terms.index,
            s = significant_terms['size']*10,
            c = significant_terms['nes'],
            cmap = 'coolwarm'
        )
        ax.set_xlabel('Enrichment score')
        ax.set_ylabel('Term')
        ax.set_title(database)
        st.pyplot(fig)

# Define the available databases
databases = [
    'GO_Biological_Process_2018',
    'GO_Cellular_Component_2018',
    'GO_Molecular_Function_2018',
    'KEGG_2019_Human',
    'KEGG_2019_Mouse',
    'KEGG_2019_Rat',
    'WikiPathways_2019_Human',
    'WikiPathways_2019_Mouse',
    'WikiPathways_2019_Rat'
]

# Create the main page of the app
st.title('Gene Set Enrichment Analysis')
st.write('Input a list of genes and select the databases to analyze.')

# Create the gene input box
genes_input = st.text_input('Genes (one per row):')

# Create the organism dropdown
organism = st.selectbox(
    'Organism:',
    ['Homo sapiens', 'Mus musculus', 'Rattus norvegicus']
)

# Create the submit button
if st.button('Submit'):
    # Process the genes input box
    genes_list = genes_input.split('\n')
    genes_list = [gene.strip() for gene in genes_list]
    genes = pd.Series(genes_list)
    # Perform GSEA and plot the results for each selected database
    for database in databases:
        perform_gsea(database, genes, organism)


