import streamlit as st
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

def perform_gsea(database, genes, organism):
    r = gp.gsea(data=genes, gene_sets=database, organism=organism, permutation_num=1000, method='signal_to_noise', processes=1, seed=1, outdir=None, no_plot=True)
    s = r.results.loc[(r.results['fdr'] < 0.05)]

    if (not s.empty):
        (f, a) = plt.subplots()
        a.scatter(s['es'], s.index, s=(s['size'] * 10), c=s['nes'], cmap='coolwarm')
        a.set_xlabel('Enrichment score')
        a.set_ylabel('Term')
        a.set_title(database)
        st.pyplot(f)

def perform_enrichr(database, genes):
    enr = gp.enrichr(gene_list=genes, gene_sets=database, background='hsapiens_gene_ensembl', outdir=None)
    enr_results = enr.results.loc[enr.results["Adjusted P-value"] < 0.05]

    if not enr_results.empty:
        ax = gp.dotplot(enr_results, column="Adjusted P-value", x='Gene_set', size=10, top_term=5, figsize=(3, 5), title=database, xticklabels_rot=45, show_ring=True, marker='o')
        st.pyplot(ax.figure)

databases = ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018', 'KEGG_2019_Human', 'KEGG_2019_Mouse', 'KEGG_2019_Rat', 'WikiPathways_2019_Human', 'WikiPathways_2019_Mouse', 'WikiPathways_2019_Rat']

st.title('Gene Set Enrichment Analysis')
st.write('Input a list of genes and select the databases to analyze.')

genes = st.text_input('Genes (one per row):')
organism = st.selectbox('Organism:', ['Homo sapiens', 'Mus musculus', 'Rattus norvegicus'])
analysis_method = st.radio("Select analysis method:", ('GSEA', 'Enrichr'))

if st.button('Submit'):
    genes = genes.split('\n')
    genes = [gene.strip() for gene in genes]
    genes = pd.Series(genes)

    for database in databases:
        if analysis_method == 'GSEA':
            perform_gsea(database, genes, organism)
        elif analysis_method == 'Enrichr':
            perform_enrichr(database, genes)
