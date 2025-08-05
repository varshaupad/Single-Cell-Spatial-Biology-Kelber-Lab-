"""
Author: Raynah Cheng and Varsha Upadhyayulla
Assignment Title: 6K RNA Scanpy Imaging and UMAPs
Date Created: 1/13/25
Date Last Modified: 1/13/25
"""
import scanpy as sc
import pandas as pd
import numpy as np

def pp(core, name):
    df = pd.read_csv(core, header = 0, low_memory = False)

    columns_removed = ['fov', 'cell_ID']

    
    df = df.drop(columns = columns_removed)
    df = df.drop(columns = [col for col in df.columns if col.startswith("Negative") or col.startswith("SystemControl")])

    standard_genes = ['EPCAM', 'CDH1', 'PTPRC', 'CD3D', 'CD3E', 'CD3E', 'CD3G', 'CD19', 'GNLY', 'NKG7', 'VCAN', 
                         'LYZ', 'ITGAM', 'CD68', 'VIM', 'FN1', 'PECAM1', 'FLT1', 'DCN', 'IL6', 'CFD', 'CD74', 'HLA-DRB', 
                         'ACTA2', 'PDGFRB', 'TAGLN', 'RGS5']


    # selected_columns = ['ACTN1', 'ACTN4', 'AKT1', 'AKT2', 'AKT3', 'ALCAM', 'ANG', 'ARHGAP5', 'BAD', 'BCAR1', 'BCL2', 
    #                   'BIRC2', 'BIRC3', 'BRAF', 'CADM1', 'CAPN2', 'CAV1', 'CAV2', 'CCND1', 'CCND2', 'CCND3', 'CD151',
    #                   'CD2', 'CD22', 'CD226', 'CD274', 'CD276', 'CD28', 'CD34', 'CD4', 'CD40', 'CD40LG', 'CD58', 'CD6', 
    #                   'CD80', 'CD86', 'CD8A', 'CD8B', 'CD99', 'CD99L2', 'CDH1', 'CDH10', 'CDH11', 'CDH12', 'CDH13', 'CDH18',
    #                   'CDH2', 'CDH3', 'CDH5', 'CDH8', 'CDH9', 'CHAD', 'CLDN1', 'CLDN3', 'CLDN4', 'CLDN5', 'CLDN7', 'CLMP', 
    #                   'CNTN1', 'CNTN2', 'COBLL1', 'COL17A1', 'COL1A1', 'COL1A2', 'COL2A1', 'COL4A1', 'COL4A2', 'COL4A3', 
    #                   'COL4A4', 'COL4A5', 'COL4A6', 'COL6A1', 'COL6A2', 'COL6A3', 'COL6A6', 'COL9A2', 'COL9A3', 'COMP',
    #                   'CRK', 'CRKL', 'CTLA4', 'CTNNA1', 'CTNNB1', 'CTNND1', 'DIAPH1', 'DOCK1', 'DPT', 'DST', 'EGF', 'EGFR', 
    #                   'ELK1', 'EMCN', 'ERBB2', 'ESAM', 'F11R', 'FLNA', 'FLNB', 'FLNC', 'FLT1', 'FN1', 'FYN', 'GLG1', 'GRB2', 
    #                   'GSK3B', 'HGF', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 
    #                   'HLA-DRA', 'HLA-E', 'HLA-F', 'HLA-G', 'HRAS', 'IBSP', 'ICAM1', 'ICAM2', 'ICAM3', 'ICOS', 'ICOSLG', 'IGF1', 
    #                   'IGF1R', 'ILK', 'ITGA1', 'ITGA2', 'ITGA2B', 'ITGA3', 'ITGA4', 'ITGA5', 'ITGA6', 'ITGA7', 'ITGA8', 'ITGA9', 
    #                   'ITGAL', 'ITGAM', 'ITGAV', 'ITGB1', 'ITGB2', 'ITGB3', 'ITGB4', 'ITGB5', 'ITGB6', 'ITGB7', 'ITGB8', 'JAM2', 
    #                   'JAM3', 'JUN', 'JUP', 'KDR', 'KRT14', 'KRT5', 'L1CAM', 'LAMA1', 'LAMA2', 'LAMA3', 'LAMA4', 'LAMA5', 'LAMB1', 
    #                 'LAMB3', 'LAMB4', 'LAMC1', 'LAMC2', 'LAMC3', 'LRRC4B', 'LRRC4C', 'MADCAM1', 'MAG', 'MAP2K1', 'MAPK1', 
    #                   'MAPK10', 'MAPK3', 'MAPK8', 'MAPK9', 'MET', 'MPZL1', 'MYL12A', 'MYL2', 'MYL7', 'MYL9', 'MYLK3', 'NCAM1', 'NECTIN1', 
    #                   'NECTIN2', 'NECTIN3', 'NECTIN4', 'NEGR1', 'NEO1', 'NFASC', 'NLGN3', 'NLGN4X', 'NLGN4Y', 'NRCAM', 'NRXN1', 'NRXN2', 
    #                   'NRXN3', 'NTNG1', 'PAK1', 'PAK2', 'PAK3', 'PAK5', 'PARD3', 'PARD6A', 'PARVB', 'PARVG', 'PDCD1', 'PDCD1LG2', 'PDGFA',
    #                    'PDGFB', 'PDGFC', 'PDGFD', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PECAM1', 'PGF', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3R1',
    #                    'PIK3R2', 'PIK3R3', 'PIP5K1C', 'PLEC', 'PODXL', 'PODXL2', 'PPP1CA', 'PPP1CC', 'PPP1R12A', 'PPP1R12B', 'PPP1R12C',
    #                    'PRKCA', 'PRKCB', 'PRKCG', 'PRKCI', 'PTEN', 'PTK2', 'PTPRC', 'PTPRF', 'PTPRM', 'PVR', 'PXN', 'RAC1', 'RAC2', 
    #                    'RAC3', 'RAF1', 'RAP1A', 'RAP1B', 'RAPGEF1', 'RASGRF1', 'RELN', 'RHOA', 'ROCK1', 'ROCK2', 'RSU1', 'SDC1', 'SDC2', 
    #                    'SDK1', 'SELE', 'SELL', 'SELP', 'SELPLG', 'SHC1', 'SHC2', 'SHC3', 'SHC4', 'SIGLEC1', 'SOS1', 'SOS2', 'SPN', 'SPP1', 
    #                    'SRC', 'THBS1', 'THBS2', 'THBS4', 'TIGIT', 'TLN1', 'TNC', 'TNN', 'TNR', 'VASP', 'VAV1', 'VAV3', 'VCAM1', 'VCAN', 
    #                    'VCL', 'VEGFA', 'VEGFB', 'VEGFC', 'VEGFD', 'VSIR', 'VTCN1', 'VTN', 'VWF', 'XIAP', 'ZYX', 'RAI14', 'PEAK1']

    #df = df[standard_genes + selected_columns]

    adata = sc.AnnData(df)
    adata.obs.index = adata.obs.index.astype(str)

    adata.var_names_make_unique()

    sc.pp.calculate_qc_metrics(adata, percent_top = None, log1p = False, inplace = True)

    sc.pp.normalize_total(adata, target_sum = 10000)
    sc.pp.log1p(adata)

    adata.raw = adata

    sc.pp.scale(adata, max_value = 10)

    sc.pp.pca(adata, svd_solver = 'arpack')
    sc.pp.neighbors(adata, n_neighbors = 99, n_pcs = 10)

    sc.tl.umap(adata)

    sc.tl.leiden(adata, key_added = 'leiden_res_0.5', resolution = 0.5, flavor = 'igraph', n_iterations = 2, random_state = 0)

    
    #sc.pl.umap(adata, color = 'leiden_res_0.5', legend_loc = "on data", save = f'_{name}.png')

    #sc.pl.umap(adata, color = genes_to_visualize, cmap = 'viridis', vmin = 0, vmax = 'p95', save = f'_{name}_markergenes.png')


    sc.tl.rank_genes_groups(adata, 'leiden_res_0.5', method = 'wilcoxon')

    

    sc.pl.rank_genes_groups(adata, n_genes = 30, sharey = False, save = f'_{name}_top30genes.png')
    #sc.pl.rank_genes_groups_heatmap(adata,  n_genes = 20, show_gene_labels = True, save = f'_{name}_top20genes.png')


    return adata




def classify_expression(adata, gene):
    median = np.median(adata[:, gene].X)
    return (adata[:, gene].X > median).flatten()

def concatenate_cell_populations(data, populations):
    concatenated_mask = np.zeros(data.shape[0], dtype=bool)
    for population in populations:
        population_mask = np.ones(data.shape[0], dtype=bool)
        for gene in marker_genes[population]:
            if gene.endswith('_pos'):
                base_gene = gene[:-4] 
                population_mask &= classify_expression(data, base_gene)
            elif gene.endswith('_neg'):
                base_gene = gene[:-4]
                population_mask &= ~classify_expression(data, base_gene)
        concatenated_mask |= population_mask
    return concatenated_mask

def create_dotplots(data, name, all_cell_types, genes_to_visualize):

    #concatenating cell definitions together
    data.obs['T cells'] = concatenate_cell_populations(data, ['T cells 1', 'T cells 2', 'T cells 3'])
    data.obs['B cells'] = concatenate_cell_populations(data, ['B cells 1', 'B cells 2', 'B cells 3'])
    data.obs['NK cells'] = concatenate_cell_populations(data, ['NK cells 1', 'NK cells 2', 'NK cells 3', 'NK cells 4', 'NK cells 5', 'NK cells 6'])
    data.obs['Monocytes'] = concatenate_cell_populations(data, ['Monocytes'])
    data.obs['Macrophages'] = concatenate_cell_populations(data, ['Macrophages'])
    data.obs['Neutrophils'] = concatenate_cell_populations(data, ['Neutrophils'])
    data.obs['Epithelial cells'] = concatenate_cell_populations(data, ['Epithelial cells'])
    data.obs['Endothelial cells'] = concatenate_cell_populations(data, ['Endothelial cells'])
    data.obs['Hybrid EMT cells'] = concatenate_cell_populations(data, ['Hybrid EMT cells'])
    data.obs['Fibroblasts (non-CAFs)'] = concatenate_cell_populations(data, ['Fibroblasts (non-CAFs)'])
    data.obs['iCAFs'] = concatenate_cell_populations(data, ['iCAFs'])
    data.obs['apCAFs'] = concatenate_cell_populations(data, ['apCAFs'])
    data.obs['myCAFs'] = concatenate_cell_populations(data, ['myCAFs 1', 'myCAFs 2'])
    data.obs['Pericytes'] = concatenate_cell_populations(data, ['Pericytes'])

    

    #sc.pl.dotplot(data, var_names = all_cell_types, groupby = 'leiden_res_0.5', save = f'{name}.png')
    #sc.pl.dotplot(data, var_names = genes_to_visualize, groupby = 'leiden_res_0.5', save = f'{name}_markergene.png')



marker_genes = {
    "T cells 1": ['EPCAM_neg', 'CDH1_neg', 'PTPRC_pos', 'CD3D_pos'], 
    "T cells 2": ['EPCAM_neg', 'CDH1_neg', 'PTPRC_pos', 'CD3E_pos'],
    "T cells 3": ['EPCAM_neg', 'CDH1_neg', 'PTPRC_pos', 'CD3G_pos'],

    "B cells 1": ['EPCAM_neg', 'CDH1_neg', 'PTPRC_pos', 'CD19_pos', 'CD3D_neg'],
    "B cells 2": ['EPCAM_neg', 'CDH1_neg', 'PTPRC_pos', 'CD19_pos', 'CD3E_neg'],
    "B cells 3": ['EPCAM_neg', 'CDH1_neg', 'PTPRC_pos', 'CD19_pos', 'CD3G_neg'],

    "NK cells 1": ['EPCAM_neg', 'CDH1_neg', 'PTPRC_pos', 'CD3D_neg', 'GNLY_pos'],
    "NK cells 2": ['EPCAM_neg', 'CDH1_neg', 'PTPRC_pos', 'CD3E_neg', 'GNLY_pos'],
    "NK cells 3": ['EPCAM_neg', 'CDH1_neg', 'PTPRC_pos', 'CD3G_neg', 'GNLY_pos'],
    "NK cells 4": ['EPCAM_neg', 'CDH1_neg', 'PTPRC_pos', 'CD3D_neg', 'NKG7_pos'],
    "NK cells 5": ['EPCAM_neg', 'CDH1_neg', 'PTPRC_pos', 'CD3E_neg', 'NKG7_pos'],
    "NK cells 6": ['EPCAM_neg', 'CDH1_neg', 'PTPRC_pos', 'CD3G_neg', 'NKG7_pos'],

    "Monocytes": ['EPCAM_neg', 'CDH1_neg', 'PTPRC_pos', 'VCAN_pos', 'LYZ_pos'],

    "Macrophages": ['EPCAM_neg', 'CDH1_neg', 'PTPRC_pos', 'ITGAM_pos', 'CD68_pos'],

    "Neutrophils": ['EPCAM_neg', 'CDH1_neg', 'PTPRC_pos', 'ITGAM_pos', 'CD68_neg'],

    "Epithelial cells": ['PTPRC_neg', 'VIM_neg', 'FN1_neg', 'EPCAM_pos', 'CDH1_pos'],

    "Endothelial cells": ['PTPRC_neg', 'PECAM1_pos', 'FLT1_pos'],

    "Hybrid EMT cells": ['PTPRC_neg', 'VIM_pos', 'CDH1_pos'],

    "Fibroblasts (non-CAFs)": ['PTPRC_neg', 'EPCAM_neg', 'CDH1_neg', 'PECAM1_neg', 'FN1_pos', 'DCN_pos', 
                               'IL6_neg', 'CFD_neg', 'CD74_neg', 'HLA-DRB_neg', 'ACTA2_neg', 'PDGFRB_neg', 'TAGLN_neg', 'RGS5_neg'],

    "iCAFs": ['PTPRC_neg', 'EPCAM_neg', 'CDH1_neg', 'PECAM1_neg', 'FN1_pos', 'DCN_pos', 'IL6_pos', 'CFD_pos',
              'CD74_neg', 'HLA-DRB_neg', 'ACTA2_neg', 'PDGFRB_neg', 'TAGLN_neg', 'RGS5_neg'],

    "apCAFs": ['PTPRC_neg', 'EPCAM_neg', 'CDH1_neg', 'PECAM1_neg', 'FN1_pos', 'DCN_pos', 'CD74_pos', 'HLA-DRB_pos',
               'IL6_neg', 'CFD_neg', 'ACTA2_neg', 'PDGFRB_neg', 'TAGLN_neg', 'RGS5_neg'],

    "myCAFs 1": ['PTPRC_neg', 'EPCAM_neg', 'CDH1_neg', 'PECAM1_neg', 'FN1_pos', 'DCN_pos', 'ACTA2_pos', 'PDGFRB_pos',
                 'IL6_neg', 'CFD_neg', 'CD74_neg', 'HLA-DRB_neg', 'TAGLN_neg', 'RGS5_neg'],
    "myCAFs 2": ['PTPRC_neg', 'EPCAM_neg', 'CDH1_neg', 'PECAM1_neg', 'FN1_pos', 'DCN_pos', 'ACTA2_pos', 'TAGLN_pos',
                 'IL6_neg', 'CFD_neg', 'CD74_neg', 'PDGFRB_neg', 'RGS5_neg'],

    "Pericytes": ['PTPRC_neg', 'EPCAM_neg', 'CDH1_neg', 'PECAM1_neg', 'FN1_pos', 'DCN_pos', 'RGS5_pos',
                  'IL6_neg', 'CFD_neg', 'CD74_neg', 'HLA-DRB_neg', 'ACTA2_neg', 'PDGFRB_neg', 'TAGLN_neg']
}


genes_to_visualize = list(set([
    #Original genes
    "TFF1", "GATA6", "KRT17", "KRT5", "KRT6A", "KRT6B", "KRT6C", "S100A2", "EPCAM", "CDH1", 
    "VIM", "CDH2", "ZEB1", "PEAK1", "ITGA1", "RAI14", "DOCK9", "INHBA", "PLK1", "FBN1", "CD36", 
]))

#from collections import OrderedDict

#genes_to_visualize = list(OrderedDict.fromkeys([
    
    #"LCN2", "SLC22A17", "MC4R", "LRP2", "ADGRL3", "BCL11B", "C1QTNF7",
    #"CCL19", "HOXB3", "HSD11B1", "LTF", "PAM", "PTPRU", "SP100",
    #"TFF2", "TIMP3", "VPREB3", "ZFHX3"

#]))





all_cell_types = ['T cells', 'B cells', 
                  'NK cells', 'Monocytes', 'Macrophages', 'Neutrophils', 
                  'Epithelial cells', 'Endothelial cells', 'Hybrid EMT cells', 
                  'Fibroblasts (non-CAFs)', 'iCAFs', 'apCAFs', 'myCAFs', 'Pericytes']

#BREAST

#HER2 Group
#her2 = pp("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_HER2.csv", "6K_RNA_HER2")
#create_dotplots(her2, "6K_RNA_HER2", all_cell_types, genes_to_visualize)

#her2_primary = pp("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_HER2_Primary.csv", "6K_RNA_HER2_Primary")
#create_dotplots(her2_primary, "6K_RNA_HER2_Primary2", all_cell_types, genes_to_visualize)

#her2_met = pp("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_HER2_Met.csv", "6K_RNA_HER2_Met")
#create_dotplots(her2_met, "6K_RNA_HER2_Met2", all_cell_types, genes_to_visualize)

#her2_primary_low = pp("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_HER2_Primary_Low_PTPRC.csv", "6K_RNA_HER2_Primary_LowPTPRC")
#create_dotplots(her2_primary_low, "6K_RNA_HER2_Primary_LowPTPRC", all_cell_types, genes_to_visualize)

#her2_primary_high = pp("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_HER2_Primary_High_PTPRC.csv", "6K_RNA_HER2_Primary_HighPTPRC")
#create_dotplots(her2_primary_high, "6K_RNA_HER2_Primary_HighPTPRC", all_cell_types, genes_to_visualize)

#her2_met_low = pp("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_HER2_Met_Low_PTPRC.csv", "6K_RNA_HER2_Met_LowPTPRC")
#create_dotplots(her2_met_low, "6K_RNA_HER2_Met_LowPTPRC", all_cell_types, genes_to_visualize)

#her2_met_high = pp("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_HER2_Met_High_PTPRC.csv", "6K_RNA_HER2_Met_HighPTPRC")
#create_dotplots(her2_met_high, "6K_RNA_HER2_Met_HighPTPRC", all_cell_types, genes_to_visualize)

#ERPR Group

#erpr_primary = pp("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_Over_or_Equal_to_40_Met.csv", "6K_RNA_ERPR_Primary")
#create_dotplots(erpr_primary, "6K_RNA_ERPR_Primary2", all_cell_types, genes_to_visualize)

#erpr_met = pp("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_ERPR_Met.csv", "6K_RNA_ERPR_Met")
#create_dotplots(erpr_met, "6K_RNA_ERPR_Met2", all_cell_types, genes_to_visualize)

#TNBC Group

#tnbc_primary = pp("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_TNBC_Primary.csv", "6K_RNA_TNBC_Primary")
#create_dotplots(tnbc_primary, "6K_RNA_TNBC_Primary2", all_cell_types, genes_to_visualize)

#tnbc_met = pp("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_TNBC_Met.csv", "6K_RNA_TNBC_Met")
#create_dotplots(tnbc_met, "6K_RNA_TNBC_Met2", all_cell_types, genes_to_visualize)


#Age-Based Groups

#over40_primary = pp("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_Over_or_Equal_to_40_Primary.csv", "6K_RNA_Over_40_Primary")
#create_dotplots(over40_primary, "6K_RNA_Over_40_Primary", all_cell_types, genes_to_visualize)

#over40_met = pp("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_Over_or_Equal_to_40_Met.csv", "6K_RNA_Over_40_Met")
#create_dotplots(over40_met, "6K_RNA_Over_40_Met", all_cell_types, genes_to_visualize)

#under40_primary = pp("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_Under_40_Primary.csv", "6K_RNA_Under_40_Primary")
#create_dotplots(under40_primary, "6K_RNA_Under_40_Primary", all_cell_types, genes_to_visualize)

#under40_met = pp("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_Under_40_Met.csv", "6K_RNA_Under_40_Met")
#create_dotplots(under40_met, "6K_RNA_Under_40_Met", all_cell_types, genes_to_visualize)



#erpr_met.obs.to_csv("/home/upadhyayullav/erprmet_obs.csv")






#Ungrouped - Over and Under 40 Primary & Met

#Uncomment Over 40's and Under 40's
#combined = over40_primary.concatenate(over40_met, under40_primary, under40_met, batch_key='sample_group')
#sc.pp.pca(combined, svd_solver='arpack')
#sc.pp.neighbors(combined, n_neighbors=99, n_pcs=10)
#sc.tl.umap(combined)
#sc.tl.leiden(combined, key_added='leiden_res_0.5', resolution=0.5, flavor='igraph', n_iterations=2)
#sc.pl.umap(combined, color='leiden_res_0.5', legend_loc="on data", save='_6K_RNA_Combined_Ungrouped2nd.png')
#sc.pl.umap(combined, color=genes_to_visualize, cmap='viridis', vmin=0, vmax='p95', save='_6K_RNA_Combined_Ungrouped2nd_markergenes.png')
#sc.tl.rank_genes_groups(combined, 'leiden_res_0.5', method='wilcoxon')
#sc.pl.rank_genes_groups(com#bined, n_genes=20, sharey=False, save='_6K_RNA_Combined_Ungrouped2nd_top20genes.png')
#sc.pl.rank_genes_groups_heatmap(combined, n_genes=20, show_gene_labels=True, save='_6K_RNA_Combined_Ungrouped2nd_top20genes.png')
#create_dotplots(combined, "6K_RNA_Combined_Ungrouped2nd", all_cell_types, genes_to_visualize)
#print(over40_primary)


#combined.uns.to_csv("/home/upadhyayullav/figures/figures/squidpy_annotated_data_var.csv")


#Primary All
#Uncomment primary age based groups
#combined_primary = over40_primary.concatenate(under40_primary, batch_key='sample_group')
#sc.pp.pca(combined_primary, svd_solver='arpack')
#sc.pp.neighbors(combined_primary, n_neighbors=99, n_pcs=10)
#sc.tl.umap(combined_primary)
#sc.tl.leiden(combined_primary, key_added='leiden_res_0.5', resolution=0.5, flavor='igraph', n_iterations=2)
#sc.pl.umap(combined_primary, color='leiden_res_0.5', legend_loc="on data", save='_Primary_All_UMAP.png')
#sc.pl.umap(combined_primary, color=genes_to_visualize, cmap='viridis', vmin=0, vmax='p95', save='_Primary_All_MarkerGenes.png')
#sc.tl.rank_genes_groups(combined_primary, 'leiden_res_0.5', method='wilcoxon')
#sc.pl.rank_genes_groups(combined_primary, n_genes=30, sharey=False, save='_Primary_All_Top30Genes.png')
#sc.pl.rank_genes_groups_heatmap(combined_primary, n_genes=30, show_gene_labels=True, save='_Primary_All_Top30Genes.png')
#create_dotplots(combined_primary, "Primary_All", all_cell_types, genes_to_visualize)


#Met All
#Uncomment met age based groups
#combined_met = over40_met.concatenate(under40_met, batch_key='sample_group')
#sc.pp.pca(combined_met, svd_solver='arpack')
#sc.pp.neighbors(combined_met, n_neighbors=99, n_pcs=10)
#sc.tl.umap(combined_met)
#sc.tl.leiden(combined_met, key_added='leiden_res_0.5', resolution=0.5, flavor='igraph', n_iterations=2)
#sc.pl.umap(combined_met, color='leiden_res_0.5', legend_loc="on data", save='_Met_All_UMAP.png')
#sc.pl.umap(combined_met, color=genes_to_visualize, cmap='viridis', vmin=0, vmax='p95', save='_Met_All_MarkerGenes.png')
#sc.tl.rank_genes_groups(combined_met, 'leiden_res_0.5', method='wilcoxon')
#sc.pl.rank_genes_groups(combined_met, n_genes=30, sharey=False, save='_Met_All_Top30Genes.png')
#sc.pl.rank_genes_groups_heatmap(combined_met, n_genes=30, show_gene_labels=True, save='_Met_All_Top30Genes.png')
#create_dotplots(combined_met, "Met_All", all_cell_types, genes_to_visualize)








#PANCREAS
#normal_nat = pp("/data/kelberj_shared/data/PA722_151_6KPanel/grouped/6K_RNA_Normal_NAT.csv", '6K_RNA_Normal_NAT')
# create_dotplots(normal_nat, '6K_RNA_Normal_NAT', all_cell_types, genes_to_visualize)

#malignant = pp("/data/kelberj_shared/data/PA722_151_6KPanel/grouped/6K_RNA_Malignant.csv", '6K_RNA_Malignant')
# create_dotplots(malignant, '6K_RNA_Malignant', all_cell_types, genes_to_visualize)

#basal_malignant = pp("/data/kelberj_shared/data/PA722_151_6KPanel/grouped/6K_RNA_Basal_Malignant.csv", '6K_RNA_Basal_Malignant')
# create_dotplots(basal_malignant, '6K_RNA_Basal_Malignant', all_cell_types, genes_to_visualize)

#classical_malignant = pp('/data/kelberj_shared/data/PA722_151_6KPanel/grouped/6K_RNA_Classical_Malignant.csv', '6K_RNA_Classical_Malignant')
# create_dotplots(classical_malignant, '6K_RNA_Classical_Malignant', all_cell_types, genes_to_visualize)

# normal = pp("/data/kelberj_shared/data/PA722_151_6KPanel/grouped/6K_RNA_Normal.csv", '6K_RNA_Normal')
# create_dotplots(normal, '6K_RNA_Normal', all_cell_types, genes_to_visualize)

# nat = pp("data/kelberj_shared/data/PA722_151_6KPanel/grouped/6K_RNA_NAT.csv", '6K_RNA_NAT')
# create_dotplots(nat, "6K_RNA_NAT", all_cell_types, genes_to_visualize)

# grade2 = pp("data/kelberj_shared/data/PA722_151_6KPanel/grouped/6K_RNA_Grade_2.csv", '6K_RNA_Grade_2')
# create_dotplots(grade2, '6K_RNA_Grade_2', all_cell_types, genes_to_visualize)

# grade3 = pp("data/kelberj_shared/data/PA722_151_6KPanel/grouped/6K_RNA_Grade_3.csv", '6K_RNA_Grade_3')
# create_dotplots(grade3, '6K_RNA_Grade_3', all_cell_types, genes_to_visualize)

# basal = pp("data/kelberj_shared/data/PA722_151_6KPanel/grouped/6K_RNA_Basal.csv", '6K_RNA_Basal')
# create_dotplots(basal, '6K_RNA_Basal', all_cell_types, genes_to_visualize)

# classical = pp("data/kelberj_shared/data/PA722_151_6KPanel/grouped/6K_RNA_Classical.csv", '6K_RNA_Classical')
# create_dotplots(classical, '6K_RNA_Classical', all_cell_types, genes_to_visualize)

over60_malignant = pp('/data/kelberj_shared/data/PA722_151_6KPanel/grouped/6K_RNA_Over60_Malignant.csv', '6K_RNA_Over60_Malignant')
#over60_malignant_clusters = over60_malignant.obs['leiden_res_0.5'].unique().tolist()

#for cluster in over60_malignant_clusters:
    #over60_malignant_df = sc.get.rank_genes_groups_df(over60_malignant, group=cluster)
    #over60_malignant_df.to_csv(f"over60_malignant_ranked_genes_cluster_{cluster}.csv", index=False)
#create_dotplots(over60_malignant, '6K_RNA_Over60_Malignant', all_cell_types, genes_to_visualize)

under60_malignant = pp('/data/kelberj_shared/data/PA722_151_6KPanel/grouped/6K_RNA_Under60_Malignant.csv', '6K_RNA_Under60_Malignant')
#under60_malignant_clusters = under60_malignant.obs['leiden_res_0.5'].unique().tolist()

#for cluster in under60_malignant_clusters:
    #under60_malignant_df = sc.get.rank_genes_groups_df(under60_malignant, group=cluster)
    #under60_malignant_df.to_csv(f"under60_malignant_ranked_genes_cluster_{cluster}.csv", index=False)
#create_dotplots(under60_malignant, '6K_RNA_Under60_Malignant', all_cell_types, genes_to_visualize)
