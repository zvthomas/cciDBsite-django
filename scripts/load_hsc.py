from tkinter import E
from interactions.models import Pathway, Receptor, Ligand, cellClas, pathwayAndCelltype, pathwayCorrelations
import pandas as pd
import numpy as np
# import all the code to properly extract relevant data from cellchat outputs


specificity_index = pd.read_csv("C:/Users/zthomas/Documents/cciDB/cciDB/cciDBsite/staticfiles/all_pathways/pathway_scale_factors.csv")

non_barcodes = ["Chondrocytes","EC-Arteriar","EC-Arteriolar","EC-Sinusoidal",
                        "Fibroblasts","MSPC-Adipo","MSPC-Osteo","Myofibroblasts","Osteo",
                        "Osteoblasts","Pericytes","Schwann-cells","Smooth-muscle", 
                       'CLP', 'lin-', 'MEP', 'CMP', 'GMP', 'MPP', 'HSC', 'B cell', 'Dendritic cells',
                        'Eo/Baso prog.', 'Ery prog.', 'Ery/Mk prog.', 'Erythroblasts', 'Gran/Mono prog.',
                        'LMPPs', 'Mk prog.', 'Mono prog.', 'Monocytes', 'NK cells', 'Neutro prog.', 'Neutrophils',
                        'T cells', 'large pre-B.', 'pro-B', 'small pre-B.']

file_names = ["Chondrocytes", "ECArteriar", "ECArteriolar", "ECSinusoidal","Fibroblasts",
                  "MSPCAdipo", "MSPCOsteo", "Myofibroblasts", "Osteo", "Osteoblasts",
                  "Pericytes", "Schwanncells", "Smoothmuscle","CLP","lin","MEP","CMP","GMP","MPP",
                  "HSC","B cell","Dendritic cells","EoBaso prog","Ery prog","EryMk prog",
                  "Erythroblasts","GranMono prog","LMPPs","Mk prog","Mono prog","Monocytes",
                  "NK cells","Neutro prog","Neutrophils","T cells","large preB","proB","small preB"]

def process_csv(filename):
    pp = pd.read_csv(filename, header  = 0)
    pp = pp.set_index('Unnamed: 0')
    pp.index.names = [None]

    pp = pp.rename(columns = {'MPP-':'MPP'}, index = {'MPP-':'MPP'})

    return pp

def process_path_file(filename, paths, thresh, celltype, df, barcodes = None):
    
    samp = process_csv(filename)
    samp = samp.T
    samp = samp.drop(columns = non_barcodes)
    if barcodes:
        samp = samp[barcodes]
    samp = samp.T
    samp = samp[[p for p in paths if p in samp.columns]]
    
    # Drop Pathways with less than 15% active HSC
    for col in samp.columns:
        zc = (samp[col] == 0).sum()
        if zc/len(samp) > ( 1 - thresh ):
            samp = samp.drop([col], axis = 1)
    
    d_to_Series = {}
    
    # Calculate Non-Zero Fraction and Non-Zero Average
    for col in samp.columns:
        zc = (samp[col] == 0).sum()
        nzf = (1 - zc/len(samp))
        nza = samp[col].mean()*len(samp)/(len(samp)-zc)
        d_to_Series[col] = [nzf, nza]
    
    if celltype == "Osteo":
        celltype = "Osteoclasts"

    ser = pd.Series(data = d_to_Series, index = d_to_Series.keys(), name = celltype)
    
    df = df.append(ser)
    
    return df

def get_cell_code(rowname):
    if rowname in ["Chondrocytes","EC-Arteriar","EC-Arteriolar","EC-Sinusoidal",
                        "Fibroblasts","MSPC-Adipo","MSPC-Osteo","Myofibroblasts","Osteoclasts",
                        "Osteoblasts","Pericytes","Schwann-cells","Smooth-muscle"] :
        return 'n'
    if rowname in ['B cell', 'Dendritic cells',
                'Erythroblasts', 'Monocytes',
                'NK cells', 'Neutrophils', 'T cells', 'large pre-B.', 'pro-B', 'small pre-B.']:
        return 'i'
    return 'p'

def drop_duplicate_rows(df):
    df = (df.reset_index().drop_duplicates(subset='index',keep='last').set_index('index').sort_index())
    return df

signal_code = {'Cell-Cell Contact':'c','Secreted Signaling':'s',"ECM-Receptor":'e'}


def run():

    # Set up pathways and cell types and their pairs

    Pathway.objects.all().delete()
    cellClas.objects.all().delete()
    pathwayAndCelltype.objects.all().delete()
    Ligand.objects.all().delete()
    Receptor.objects.all().delete()
    pathwayCorrelations.objects.all().delete()

    path_annot = pd.read_csv('P:/zthomas/Intercellular Interactions Project/cluster_analysis/path_annotation.csv', sep = ',')
    path_annot = path_annot.set_index("paths")
    path_annot.index.names = [None]
    


    for HSPC_TYPE in ['hsc', 'clp', 'cmp', 'gmp', 'mep', 'mpp']:

        if HSPC_TYPE == 'hsc':
            path_path = "P:/zthomas/Intercellular Interactions Project/cluster_analysis/general_niche_hscMetacell/all_pathways/"
        else:
            path_path = "P:/zthomas/Intercellular Interactions Project/cluster_analysis/general_niche_" + HSPC_TYPE + "Metacell/all_pathways/"

        filepathS_names = [path_path + 'singleCell_to_' + f + ".csv" for f in file_names]
        filepathR_names = [path_path + f + "_to_singleCell.csv" for f in file_names]
        signal_code = {'Cell-Cell Contact':'c','Secreted Signaling':'s',"ECM-Receptor":'e'}
        
        for signal_type in ['Cell-Cell Contact','Secreted Signaling',"ECM-Receptor"]:
            
            paths = list(path_annot[path_annot['annot'] == signal_type].index)

            NZFA_S = pd.DataFrame()
            NZFA_R = pd.DataFrame()

            for i in range(len(filepathS_names)):
                NZFA_S = process_path_file(filepathS_names[i], paths, 0, non_barcodes[i], NZFA_S)
            NZFA_S = drop_duplicate_rows(NZFA_S)

            for i in range(len(filepathR_names)):
                NZFA_R = process_path_file(filepathR_names[i], paths, 0, non_barcodes[i], NZFA_R)
            NZFA_R = drop_duplicate_rows(NZFA_R)

            for row in NZFA_S.index:
                for col in NZFA_S.columns:
                    nzfa = NZFA_S.loc[row, col]
                    
                    if nzfa[0] == 0.0:
                        nzfa[1] = 0.0
                    
                    si = specificity_index[specificity_index['Pathway'] == col]
                    si = si[si['Direction'] == 'sending']
                    if len(si) > 0:
                        si_s = float(list(si['max(LR) scale factor'])[0])
                    else:
                        si_s = 0
                    
                    si = specificity_index[specificity_index['Pathway'] == col]
                    si = si[si['Direction'] == 'receiving']
                    if len(si) > 0:
                        si_r = float(list(si['max(LR) scale factor'])[0])
                    else:
                        si_r = 0
                    
                    
                    
                    pathway, _ = Pathway.objects.get_or_create(name = col,
                                                            interaction_type = signal_code[signal_type],
                                                            specificity_index_s = si_s,
                                                            specificity_index_r = si_r,
                                                            evidences = "")
                    
                    ct, _ = cellClas.objects.get_or_create(name = row, cell_type=get_cell_code(row))
                    ct.sendingPathways.add(pathway)


                    

                    pathwayAndCelltype.objects.get_or_create(pathway = pathway,
                                                    celltype = ct,
                                                    hscPercent = nzfa[0],
                                                    averageScore = nzfa[1],
                                                    specificity_index = si_s,
                                                    sorr = 's',
                                                    hspc_type = HSPC_TYPE)

            
            for row in NZFA_R.index:
                for col in NZFA_R.columns:
                    nzfa = NZFA_R.loc[row, col]
                    
                    if nzfa[0] == 0.0:
                        nzfa[1] = 0.0


                    si = specificity_index[specificity_index['Pathway'] == col]
                    si = si[si['Direction'] == 'sending']
                    if len(si) > 0:
                        si_s = float(list(si['max(LR) scale factor'])[0])
                    else:
                        si_s = 0
                    
                    si = specificity_index[specificity_index['Pathway'] == col]
                    si = si[si['Direction'] == 'receiving']
                    if len(si) > 0:
                        si_r = float(list(si['max(LR) scale factor'])[0])
                    else:
                        si_r = 0


                    pathway, _ = Pathway.objects.get_or_create(name = col,
                                                            interaction_type = signal_code[signal_type],
                                                            specificity_index_s = si_s,
                                                            specificity_index_r = si_r,
                                                            evidences = "")
                    
                    ct, _ = cellClas.objects.get_or_create(name = row, cell_type=get_cell_code(row))
                    ct.receivingPathways.add(pathway)


                    pathwayAndCelltype.objects.get_or_create(pathway = pathway,
                                                    celltype = ct,
                                                    hscPercent = nzfa[0],
                                                    averageScore = nzfa[1],
                                                    specificity_index = si_r,
                                                    sorr = 't',
                                                    hspc_type = HSPC_TYPE)
                
    
    
    
    
    # Set up ligands and receptors

    
    interaction_genes = process_csv("P:/zthomas/Intercellular Interactions Project/interaction_genes.csv")
    
    for row in interaction_genes.index:

        pairs = interaction_genes.loc[row, ]
        evidence = pairs['evidence']
        ligand, _ = Ligand.objects.get_or_create(name = pairs['ligand'])
        receptor, _ = Receptor.objects.get_or_create(name = pairs['receptor'])
        pathway, _ = Pathway.objects.get_or_create(name = pairs['pathway_name'], interaction_type = signal_code[pairs['annotation']])
        ligand.receptors.add(receptor)
        ligand.pathways.add(pathway)

        receptor.ligands.add(ligand)
        receptor.pathways.add(pathway)

        pathway.ligands.add(ligand)
        pathway.receptors.add(receptor)

        pathway.evidences = pathway.evidences + evidence + ";"
        pathway.save()

        receptor.evidences = receptor.evidences + evidence + ";"
        ligand.evidences = ligand.evidences + evidence + ";"
        receptor.save()
        ligand.save()

    # Set up z-score correlations
    
    for HSPC_TYPE in ['hsc', 'clp', 'cmp', 'gmp', 'mep', 'mpp']:

        pos_correlations = pd.read_csv("P:/zthomas/Intercellular Interactions Project/cluster_analysis/pathway_correlations/" + HSPC_TYPE + "Metacell_meanpos_corr_zscore_500_v1.csv")
        neg_correlations = pd.read_csv("P:/zthomas/Intercellular Interactions Project/cluster_analysis/pathway_correlations/" + HSPC_TYPE + "Metacell_meanneg_corr_zscore_500_v1.csv")
        
        # change send and receive to source and target
        name_change = {'s':'s', 'r':'t'}
        for i in [pos_correlations, neg_correlations]:
                for row in i.index:
                    cor = i.loc[row,]

                    p1 = cor['Pathway 1'].split()[1:][0]
                    p1a = name_change[cor['Pathway 1'].split()[0][0]]
                    p2 = cor['Pathway 2'].split()[1:][0]
                    p2a = name_change[cor['Pathway 2'].split()[0][0]]

                    corr = cor['Correlation']
                    pval = cor['Adj. P-val']
                    
                    pathway1, _ = Pathway.objects.get_or_create(name = p1)
                    pathway2, _ = Pathway.objects.get_or_create(name = p2)

                    pathwayCorrelations.objects.get_or_create(pathway1 = pathway1,
                                                            pathway2 = pathway2,
                                                            p1a = p1a,
                                                            p2a = p2a,
                                                            correlation = corr,
                                                            pval = pval,
                                                            hspc_type = HSPC_TYPE)

    