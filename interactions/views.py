from django.shortcuts import render
from plotly.offline import plot
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from django.shortcuts import get_object_or_404
import plotly.express as px
import networkx as nx
import json
import pickle
import numpy as np
from matplotlib import cm
import matplotlib as mpl
from .models import Pathway, Receptor, Ligand, cellClas, pathwayAndCelltype, pathwayCorrelations
from django.views import generic
from django.contrib.staticfiles.storage import staticfiles_storage
import pandas as pd
# Create your views here.

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 11:45:07 2020
@author: zthomas
"""

greys = cm.get_cmap('Greys', 10000)

cell_type_proportions = {'CLP': 0.074845578,
 'MPP': 0.391143494,
 'HSC':  0.00403241,
 'CMP':  0.779191511,
 'GMP':  0.384841643,
 'MEP':  0.148679469,
 'B cell': 8.081785369,
 'Monocytes': 3.068175516,
 'T cells': 2.335698011,
 'Chondrocytes': 0.106393443,
 'EC-Arteriar':  0.014297933,
 'EC-Arteriolar':  0.029577097,
 'EC-Sinusoidal':  0.081722499,
 'Fibroblasts':  0.156996911,
 'MSPC-Adipo':  0.068125445,
 'MSPC-Osteo':  0.020325493,
 'Myofibroblasts': 0.002803516,
 'Osteoclasts': 0.076115467,
 'Osteoblasts':  0.026353053,
 'Pericytes':  0.005326681,
 'Schwann-cells':  0.001962461,
 'Dendritic cells':  0.47,
 'Erythroblasts': 28.036}

KLS_spatial_pvalue = {'Neutrophils': 0.966643096, 'Erythroblasts': 0.587418233, 'B cell': 1.099775133, 'EC-Arteriar': 1.87721653,
                      'EC-Sinusoidal': 3.759205485, 'MSPC-Adipo': 6.381276781, 'MSPC-Osteo': 6.381276781, 'T cells': 1.274467619,
                      'Dendritic cells': 1.134701735, 'GMP': 2.118434338, 'MEP': 2.141933006, 'Monocytes': 1.846794428,
                      'HSC': 12.83673889, 'MPP': 12.83673889, 'CMP': 1.276287817, 'CLP': 0.647591814}

MEP_spatial_pvalue = {'Neutrophils': 1.134389547, 'Erythroblasts': 0.687436463, 'B cell': 1.219862119, 'EC-Arteriar': 1.066141032,
                      'EC-Sinusoidal': 1.350647685, 'MSPC-Adipo': 1.198134108, 'MSPC-Osteo': 1.198134108, 'T cells': 1.446986139,
                      'Dendritic cells': 1.509321242, 'GMP': 1.568872768, 'MEP': 3.588684567, 'Monocytes': 1.701691206,
                      'HSC': 2.203514657, 'MPP': 2.203514657, 'CMP': 1.706437938, 'CLP': 0.869342193}

GMP_spatial_pvalue = {'Neutrophils': 0.948266884, 'Erythroblasts': 0.854901519, 'B cell': 1.040639476, 'EC-Arteriar': 1.124449318,
                      'EC-Sinusoidal': 1.439652345, 'MSPC-Adipo': 1.428992115, 'MSPC-Osteo': 1.428992115, 'T cells': 1.279855198,
                      'Dendritic cells': 1.589494537, 'GMP': 6.878198728, 'MEP': 1.571141937, 'Monocytes': 1.611571305,
                      'HSC': 2.258687428, 'MPP': 2.258687428, 'CMP': 2.045440578, 'CLP': 0.982363337}

CMP_spatial_pvalue = {'Neutrophils': 0.964586022, 'Erythroblasts': 0.80978427, 'B cell': 1.146620302, 'EC-Arteriar': 1.265054034,
                      'EC-Sinusoidal': 1.166467984, 'MSPC-Adipo': 1.045712139, 'MSPC-Osteo': 1.045712139, 'T cells': 1.168804887,
                      'Dendritic cells': 1.114616563, 'GMP': 2.458688788, 'MEP': 1.831683131, 'Monocytes': 1.2825187,
                      'HSC': 2.185727581, 'MPP': 2.185727581, 'CMP': 6.116084621, 'CLP': 1.426587295}

CLP_spatial_pvalue = {'Neutrophils': 0.990714824, 'Erythroblasts': 0.911156383, 'B cell': 1.392372321, 'EC-Arteriar': 1.052534173,
                      'EC-Sinusoidal': 0.814265403, 'MSPC-Adipo': 0.637125047, 'MSPC-Osteo': 0.637125047, 'T cells': 0.942937992,
                      'Dendritic cells': 1.435749445, 'GMP': 0.88959256, 'MEP': 0.989320166, 'Monocytes': 1.250350075,
                      'HSC': 0.935126203, 'MPP': 0.935126203, 'CMP': 0.970787627, 'CLP': 3.676524164}

KLS_spatial_iniche = {'Neutrophils': 0.870513923, 'Erythroblasts': 0.594584591, 'B cell': 1.085288232, 'EC-Arteriar': 2.221900637,
                      'EC-Sinusoidal': 5.123134068, 'MSPC-Adipo': 8.235350998,'MSPC-Osteo': 8.235350998,'T cells': 1.012093793, 'Dendritic cells': 0.912145908,
                      'GMP': 3.122276902, 'MEP': 2.547517799, 'Monocytes': 1.688492404, 'HSC': 26.79845834,'MPP': 26.79845834, 'CMP': 1.809825199, 'CLP': 0.681595783}
MEP_spatial_iniche = {'Neutrophils': 1.027017169, 'Erythroblasts': 0.623274061, 'B cell': 1.103992649, 'EC-Arteriar': 1.262264693,
                      'EC-Sinusoidal': 1.839055505, 'MSPC-Adipo': 2.335190778,'MSPC-Osteo': 2.335190778,'T cells': 1.276282533, 'Dendritic cells': 1.534602872,
                      'GMP': 5.918830164, 'MEP': 4.526600526, 'Monocytes': 1.469816001, 'HSC': 6.615499939,'MPP': 6.615499939, 'CMP': 4.724198954, 'CLP': 0.956658488}
GMP_spatial_iniche = {'Neutrophils': 0.893082276, 'Erythroblasts': 0.668647871, 'B cell': 0.992036277, 'EC-Arteriar': 1.547088954,
                      'EC-Sinusoidal': 2.719085464, 'MSPC-Adipo': 3.846213041,'MSPC-Osteo': 3.846213041,'T cells': 0.976080972, 'Dendritic cells': 1.377448191,
                      'GMP': 10.07540326, 'MEP': 2.254418104, 'Monocytes': 1.6447961, 'HSC': 10.71426174,'MPP': 10.71426174, 'CMP': 5.583027099, 'CLP': 0.854510698}
CMP_spatial_iniche = {'Neutrophils': 0.900260883, 'Erythroblasts': 0.611567155, 'B cell': 1.002615358, 'EC-Arteriar': 2.074536174,
                      'EC-Sinusoidal': 3.986705432, 'MSPC-Adipo': 6.512222761,'MSPC-Osteo': 6.512222761,'T cells': 0.960071791, 'Dendritic cells': 0.971940154,
                      'GMP': 2.957201498, 'MEP': 2.874064379, 'Monocytes': 1.554744457, 'HSC': 20.58613684,'MPP': 20.58613684, 'CMP': 9.561074757, 'CLP': 0.878619619}
CLP_spatial_iniche = {'Neutrophils': 1.022243378, 'Erythroblasts': 0.646840852, 'B cell': 1.219614677, 'EC-Arteriar': 1.230849931,
                      'EC-Sinusoidal': 1.70212692, 'MSPC-Adipo': 2.137946093,'MSPC-Osteo':2.137946093, 'T cells': 1.193537649, 'Dendritic cells': 1.528236542,
                      'GMP': 5.358391301, 'MEP': 4.145741274, 'Monocytes': 1.379459244, 'HSC': 5.966991466,'MPP': 5.966991466, 'CMP': 4.367956317, 'CLP': 2.391773299}

pheno_dicts = {'hsc':KLS_spatial_pvalue, 'mpp':KLS_spatial_pvalue,
               'mep':MEP_spatial_pvalue, 'gmp':GMP_spatial_pvalue,
               'cmp': CMP_spatial_pvalue, 'clp': CLP_spatial_pvalue}

iniche_dicts = {'hsc':KLS_spatial_iniche, 'mpp':KLS_spatial_iniche,
               'mep':MEP_spatial_iniche, 'gmp':GMP_spatial_iniche,
               'cmp': CMP_spatial_iniche, 'clp': CLP_spatial_iniche}

r_pickle_path = staticfiles_storage.path("all_pathways/receiving_pathways_dictionary_db.pkl")
s_pickle_path = staticfiles_storage.path("all_pathways/sending_pathways_dictionary_db.pkl")
r_pickle_celltype = staticfiles_storage.path("all_pathways/receiving_celltypes_dictionary_db.pkl")
s_pickle_celltype = staticfiles_storage.path("all_pathways/sending_celltypes_dictionary_db.pkl")

specificity_index = staticfiles_storage.path("all_pathways/hspc_specificity_index_v4.csv")
specificity_index = pd.read_csv(specificity_index, index_col = 0)

with open(r_pickle_path, 'rb') as handle:
    receiving_pathways_dictionary_db = pickle.load(handle)

with open(s_pickle_path, 'rb') as handle:
    sending_pathways_dictionary_db = pickle.load(handle)

with open(r_pickle_celltype, 'rb') as handle:
    receiving_celltypes_dictionary_db = pickle.load(handle)

with open(s_pickle_celltype, 'rb') as handle:
    sending_celltypes_dictionary_db = pickle.load(handle)


import math
from typing import List
from itertools import chain

# n is number of points to place evenly in a circle centered around (0,0)
def circle(n):
    x = []
    y = []
    for i in range(0,n):
        angle = (i/n)*math.pi*2
        x.append(math.cos(angle))
        y.append(math.sin(angle))
    
    return x, y

# Start and end are lists defining start and end points
# Edge x and y are lists used to construct the graph
# arrowAngle and arrowLength define properties of the arrowhead
# arrowPos is None, 'middle' or 'end' based on where on the edge you want the arrow to appear
# arrowLength is the length of the arrowhead
# arrowAngle is the angle in degrees that the arrowhead makes with the edge
# dotSize is the plotly scatter dot size you are using (used to even out line spacing when you have a mix of edge lengths)
def addEdge(start, end, edge_x, edge_y, lengthFrac=1, arrowPos = None, arrowLength=0.025, arrowAngle = 30, dotSize=20):

    # Get start and end cartesian coordinates
    x0, y0 = start
    x1, y1 = end

    # Incorporate the fraction of this segment covered by a dot into total reduction
    length = math.sqrt( (x1-x0)**2 + (y1-y0)**2 )
    dotSizeConversion = .0565/20 # length units per dot size
    convertedDotDiameter = dotSize * dotSizeConversion
    lengthFracReduction = convertedDotDiameter / length
    lengthFrac = lengthFrac - lengthFracReduction

    # If the line segment should not cover the entire distance, get actual start and end coords
    skipX = (x1-x0)*(1-lengthFrac)
    skipY = (y1-y0)*(1-lengthFrac)
    x0 = x0 + skipX/2
    x1 = x1 - skipX/2
    y0 = y0 + skipY/2
    y1 = y1 - skipY/2

    # Append line corresponding to the edge
    edge_x.append(x0)
    edge_x.append(x1)
    edge_x.append(None) # Prevents a line being drawn from end of this edge to start of next edge
    edge_y.append(y0)
    edge_y.append(y1)
    edge_y.append(None)

    # Draw arrow
    if not arrowPos == None:

        # Find the point of the arrow; assume is at end unless told middle
        pointx = x1
        pointy = y1

        eta = math.degrees(math.atan((x1-x0)/(y1-y0))) if y1!=y0 else 90.0

        if arrowPos == 'middle' or arrowPos == 'mid':
            pointx = x0 + (x1-x0)/2
            pointy = y0 + (y1-y0)/2

        # Find the directions the arrows are pointing
        signx = (x1-x0)/abs(x1-x0) if x1!=x0 else +1    #verify this once
        signy = (y1-y0)/abs(y1-y0) if y1!=y0 else +1    #verified

        # Append first arrowhead
        dx = arrowLength * math.sin(math.radians(eta + arrowAngle))
        dy = arrowLength * math.cos(math.radians(eta + arrowAngle))
        edge_x.append(pointx)
        edge_x.append(pointx - signx**2 * signy * dx)
        edge_x.append(None)
        edge_y.append(pointy)
        edge_y.append(pointy - signx**2 * signy * dy)
        edge_y.append(None)

        # And second arrowhead
        dx = arrowLength * math.sin(math.radians(eta - arrowAngle))
        dy = arrowLength * math.cos(math.radians(eta - arrowAngle))
        edge_x.append(pointx)
        edge_x.append(pointx - signx**2 * signy * dx)
        edge_x.append(None)
        edge_y.append(pointy)
        edge_y.append(pointy - signx**2 * signy * dy)
        edge_y.append(None)


    return edge_x, edge_y

def add_arrows(source_x: List[float], target_x: List[float], source_y: List[float], target_y: List[float],
               arrowLength=0.025, arrowAngle=30):
    pointx = list(map(lambda x: x[0] + (x[1] - x[0]) / 2, zip(source_x, target_x)))
    pointy = list(map(lambda x: x[0] + (x[1] - x[0]) / 2, zip(source_y, target_y)))
    etas = list(map(lambda x: math.degrees(math.atan((x[1] - x[0]) / (x[3] - x[2]))),
                    zip(source_x, target_x, source_y, target_y)))

    signx = list(map(lambda x: (x[1] - x[0]) / abs(x[1] - x[0]), zip(source_x, target_x)))
    signy = list(map(lambda x: (x[1] - x[0]) / abs(x[1] - x[0]), zip(source_y, target_y)))

    dx = list(map(lambda x: arrowLength * math.sin(math.radians(x + arrowAngle)), etas))
    dy = list(map(lambda x: arrowLength * math.cos(math.radians(x + arrowAngle)), etas))
    none_spacer = [None for _ in range(len(pointx))]
    arrow_line_x = list(map(lambda x: x[0] - x[1] ** 2 * x[2] * x[3], zip(pointx, signx, signy, dx)))
    arrow_line_y = list(map(lambda x: x[0] - x[1] ** 2 * x[2] * x[3], zip(pointy, signx, signy, dy)))

    arrow_line_1x_coords = list(chain(*zip(pointx, arrow_line_x, none_spacer)))
    arrow_line_1y_coords = list(chain(*zip(pointy, arrow_line_y, none_spacer)))

    dx = list(map(lambda x: arrowLength * math.sin(math.radians(x - arrowAngle)), etas))
    dy = list(map(lambda x: arrowLength * math.cos(math.radians(x - arrowAngle)), etas))
    none_spacer = [None for _ in range(len(pointx))]
    arrow_line_x = list(map(lambda x: x[0] - x[1] ** 2 * x[2] * x[3], zip(pointx, signx, signy, dx)))
    arrow_line_y = list(map(lambda x: x[0] - x[1] ** 2 * x[2] * x[3], zip(pointy, signx, signy, dy)))

    arrow_line_2x_coords = list(chain(*zip(pointx, arrow_line_x, none_spacer)))
    arrow_line_2y_coords = list(chain(*zip(pointy, arrow_line_y, none_spacer)))

    x_arrows = arrow_line_1x_coords + arrow_line_2x_coords
    y_arrows = arrow_line_1y_coords + arrow_line_2y_coords

    return x_arrows, y_arrows


def getHSPCcookie(request):

    hspcType = request.COOKIES.get("hspcChoice")

    if not hspcType:
        hspcType = 'hsc'

    return hspcType

def getCorrelationcookie(request):

    corrName = request.COOKIES.get("correlationName")

    return corrName

def getPos_or_Negcookie(request):

    pos_or_neg = request.COOKIES.get("pos_or_neg")


    return pos_or_neg

def getCorrectionCookie(request):

    correction = request.COOKIES.get("correctionChoice")

    if not correction:
        correction = 'none'
    
    return correction

def normalize_dot_size(min_size, max_size, val, min_p, max_p):
    val = val*(max_size-min_size)
    val = val + min_size

    return val

def get_ec(lab, ct_or_p):

    cdict = {'n':'#D12626', 'i':'#51C206', 'p':'#0855CA'}
    if ct_or_p != 'cts':
        cdict = {'s':'#056E5D', 'c':'#7905D7', 'e':'#202A25'}

    return cdict[lab]

def index(request):
    """View function for home page"""

    # Generate counts of pathways and cell types
    num_pathways = Pathway.objects.filter().count()
    num_ss = Pathway.objects.filter(interaction_type="s").count()
    num_ccc = Pathway.objects.filter(interaction_type="c").count()
    num_ecmr = Pathway.objects.filter(interaction_type="e").count()


    num_cell_types = cellClas.objects.all().count()
    num_n = cellClas.objects.filter(cell_type="n").count()
    num_i = cellClas.objects.filter(cell_type="i").count()
    num_p = cellClas.objects.filter(cell_type="p").count()

    # cell type - pathway paris
    num_celltype_pathway_pairs = pathwayAndCelltype.objects.filter(hscPercent__gt=0.049).count()
    num_ligands = Ligand.objects.count()
    num_receptors = Receptor.objects.count()

    num_correlations = pathwayCorrelations.objects.count()
    pos_correlations = pathwayCorrelations.objects.filter(correlation__gt=0).count()
    neg_correlations = pathwayCorrelations.objects.filter(correlation__lt=0).count()

    for_table = pathwayAndCelltype.objects.all().order_by('-averageScore')[:20]

    context = {
        'num_pathways' : num_pathways,
        'num_cell_types' : num_cell_types,
        'num_celltype_pathway_pairs' : num_celltype_pathway_pairs,
        'num_ss': num_ss,
        'num_ccc': num_ccc,
        'num_ecmr': num_ecmr,
        'num_n': num_n,
        'num_i': num_i,
        'num_p': num_p,
        'num_ligands': num_ligands,
        'num_receptors': num_receptors,
        'num_correlations' : num_correlations,
        'pos_correlations' : pos_correlations,
        'neg_correlations' : neg_correlations,
        'for_table' : for_table
    }

    # Render the HTML template index.html with the data in the context variable
    return render(request, 'index.html', context = context)

def help(request):
    return render(request, 'about.html')

def downloadPage(request):
    return render(request, 'download.html')

def get_evidence_list(evidences):
    
    evidenceList = evidences.split(";")
        
    evidenceList.pop()
        
    for evidence in evidenceList:
        if ", " in evidence:
            e = evidence.split(", ")
            evidenceList.remove(evidence)
            for item in e:
                evidenceList.append(item)
        
    pmids = {}
    keggs = {}
    pmcs = {}
    
    for i in evidenceList:
        i = i.replace(" ", "")
        if "PMID" in i:
           pmids[i] = "https://pubmed.ncbi.nlm.nih.gov/" + i.split(':')[1]
        elif "KEGG" in i:
            keggs[i] = "https://www.genome.jp/dbget-bin/www_bget?pathway:" + i.split(':')[-1]
        elif "PMC" in i:
            pmcs[i] = "https://www.ncbi.nlm.nih.gov/pmc/articles/" + i
    
    return pmids, keggs, pmcs

def correlations(request):

    hspcType = getHSPCcookie(request)

    ligand_list = Ligand.objects.all()
    genemap, paths = make_gene_map(ligand_list)
    correlations = pathwayCorrelations.objects.filter(hspc_type = hspcType)
    plot_div = make_net_graph_corr(correlations)
    
    # write some function that parses it like the make_net_graph_JSON function to pass it to render...
    # in the html page, write some javascript that filters the graph based on the form submission
    
    context = {'genemap':genemap, 'paths':paths, 'plot_div' : plot_div, 'hspcType' : hspcType.upper()}
    
    return render(request, 'correlations.html', context)

def process_csv(filename):
    pp = pd.read_csv(filename, header  = 0)
    pp = pp.set_index('Unnamed: 0')
    pp.index.names = [None]
    return pp

def get_specific_interaction(parent_path, cell_type, sorr):
    
    nbc = ["Chondrocytes","EC-Arteriar","EC-Arteriolar","EC-Sinusoidal",
                        "Fibroblasts","MSPC-Adipo","MSPC-Osteo","Myofibroblasts","Osteo",
                        "Osteoblasts","Pericytes","Schwann-cells","Smooth-muscle", 
                       'CLP', 'lin-', 'MEP', 'CMP', 'GMP', 'MPP', 'HSC', 'B cell', 'Dendritic cells',
                        'Eo/Baso prog.', 'Ery prog.', 'Ery/Mk prog.', 'Erythroblasts', 'Gran/Mono prog.',
                        'LMPPs', 'Mk prog.', 'Mono prog.', 'Monocytes', 'NK cells', 'Neutro prog.', 'Neutrophils',
                        'T cells', 'large pre-B.', 'pro-B', 'small pre-B.']
    
    if sorr == 's':
        filename = 'singleCell_to_' + cell_type
    else:
        filename = cell_type + "_to_singleCell"
    
    interaction_df = process_csv(parent_path + "/" +  filename + '.csv')
    
    interaction_df = interaction_df.T
    interaction_df = interaction_df.drop(columns=nbc)
    
    return interaction_df.T

def pathway_vector(hspc, pathway, direction):
        
    parent_path = staticfiles_storage.path("all_pathways/all_pathways_" + hspc + "/")
    
    heat_vector = []
    
    
    file_names = [
                  "Chondrocytes", "ECArteriar", "ECArteriolar", "ECSinusoidal","Fibroblasts",
                  "MSPCAdipo", "MSPCOsteo", "Myofibroblasts", "Osteo", "Osteoblasts",
                  "Pericytes", "Schwanncells", "Smoothmuscle",
                  "CLP","lin","MEP","CMP","GMP","MPP","HSC",
                  "B cell","Dendritic cells",
                  "EoBaso prog","Ery prog","EryMk prog","Erythroblasts","GranMono prog",
                  "LMPPs","Mk prog","Mono prog","Monocytes","NK cells","Neutro prog",
                  
                  "Neutrophils","T cells","large preB","proB","small preB"
                 ]
    
    for i in file_names:
        t = get_specific_interaction(parent_path, i, direction)
        if pathway in t.columns:
            heat_vector.append(t[pathway].sum())
        else:
            heat_vector.append(0)
    return heat_vector



def hspc_heatmap(pathway, direction, hspcs = ['hsc', 'mpp', 'clp', 'cmp', 'gmp', 'mep']):
    
    file_names = [
                  "Chondrocytes", "ECArteriar", "ECArteriolar", "ECSinusoidal","Fibroblasts",
                  "MSPCAdipo", "MSPCOsteo", "Myofibroblasts", "Osteo", "Osteoblasts",
                  "Pericytes", "Schwanncells", "Smoothmuscle",
                  "CLP","lin","MEP","CMP","GMP","MPP","HSC",
                  "B cell","Dendritic cells",
                  "EoBaso prog","Ery prog","EryMk prog","Erythroblasts","GranMono prog",
                  "LMPPs","Mk prog","Mono prog","Monocytes","NK cells","Neutro prog",
                  
                  "Neutrophils","T cells","large preB","proB","small preB"
                 ]
    
    heat_mat = []
    
    for hspc in hspcs:
        heat_mat.append(pathway_vector(hspc, pathway, direction))
    
    df = pd.DataFrame(heat_mat,columns = file_names)
    df.index = ['HSC', 'MPP', 'CLP', 'CMP', 'GMP', 'MEP']
    df = df.T
    
    df = df/(df.max().max())
    
    return df


def correlations_heatmap(request):

    hspcType = getHSPCcookie(request)

    guess = staticfiles_storage.path('download_files/' + hspcType + '_correlations_pivot.csv')
    correlation_pivot_both = pd.read_csv(guess, index_col = 0)

    correlation_pivot_both_values = correlation_pivot_both.values.tolist()
    for i in range(len(correlation_pivot_both_values)):
        correlation_pivot_both_values[i][i] = None

    # <a href='www.google.com'>heheLINK</a>

    # Create the plot
    fig = go.Figure(data = go.Heatmap(z = correlation_pivot_both_values[::-1],
                                      x = correlation_pivot_both.columns,
                                      y = correlation_pivot_both.columns[::-1],
                                      colorscale = 'RdBu_r',
                                      zmid = 0,
                                      xgap = 0.5,
                                      ygap = 0.5,
                                      hovertemplate = "%{x} %{y}<br>Spearman Rank Corr: %{z}<extra></extra>"))
    fig['layout']['xaxis']['scaleanchor']='y'
    fig.update_layout(scene = go.layout.Scene(aspectratio = {'x':1, 'y':1}),
                      height = 800,
                      plot_bgcolor='rgba(0,0,0,0)')

    # Embed the plot in an HTML div tag
    heatmap_div = plot(fig, output_type="div",)

    max_coords = correlation_pivot_both.stack().idxmax()
    min_coords = correlation_pivot_both.stack().idxmin()

    max_sentence = max_coords[0] + " and " + max_coords[1]
    min_sentence = min_coords[0] + " and " + min_coords[1]

    context: dict = {'bar_plot': heatmap_div,
                     'max_sentence':max_sentence,
                     'min_sentence':min_sentence}

    return render(request, 'correlation_heatmap.html', context)

def signalingNetworkPage(request):

    return render(request, 'signalingNetworks.html')


def plot_correlation_scatter(request, pos_neg):
    
    hspcType = getHSPCcookie(request)

    correlation_name = request.GET.get('correlation_name')

    title = correlation_name
    if 'Source' in correlation_name:
        correlation_name = correlation_name.replace("Source", 'send')
    if 'Target' in correlation_name:
        correlation_name = correlation_name.replace("Target", 'receive')

    guess = staticfiles_storage.path("download_files/all_"+pos_neg+"_correlations_" + hspcType + "Metacell.csv")
    basis = pd.read_csv(guess, index_col = 0)

    df_max_scaled = staticfiles_storage.path('download_files/' + hspcType + 'Metacell_interactions_scaled.csv')
    df_max_scaled = pd.read_csv(df_max_scaled, index_col = 0)
    
    specific_corr = basis[basis['TrueName'] == correlation_name]
    if len(specific_corr) == 0:
        correlation_name_split = correlation_name.split(" ")
        correlation_name_split_part1 = " ".join(correlation_name_split[:2])
        correlation_name_split_part2 = " ".join(correlation_name_split[2:])
        
        correlation_name = correlation_name_split_part2 + " " + correlation_name_split_part1
        specific_corr = basis[basis['TrueName'] == correlation_name]

    
    pathway1 = pd.DataFrame()
    pathway2 = pd.DataFrame()    
    
    pathwaySplit = correlation_name.split(" ")
    pathway_name1 = " " + pathwaySplit[1]
    pathway_name2 = " " + pathwaySplit[3]
    
    for i in specific_corr['Pathways']:
        pathway_comps = i.split(" and ")
        
        if pathway_name1 in pathway_comps[0]:
            pathway1[pathway_comps[0]] = df_max_scaled[pathway_comps[0]]
            pathway2[pathway_comps[1]] = df_max_scaled[pathway_comps[1]]
        else:
            pathway1[pathway_comps[1]] = df_max_scaled[pathway_comps[1]]
            pathway2[pathway_comps[0]] = df_max_scaled[pathway_comps[0]]
    
    pathway_combo1 = pathway1.mean(axis = 1) 
    pathway_combo2 = pathway2.mean(axis = 1) 
    
    toPlot = pd.DataFrame()

    if pathwaySplit[0] == 'send':
        x_axis = pathway_name1 + ' source'
    else:
        x_axis = pathway_name1 + ' target'

    if pathwaySplit[2] == 'send':
        y_axis = pathway_name2 + ' source'
    else:
        y_axis = pathway_name2 + ' target'


    toPlot[x_axis] = pathway_combo1
    toPlot[y_axis] = pathway_combo2
    
    linecolor = 'blue'
    if pos_neg == 'positive':
        linecolor = 'red'
    
    fig = go.Figure()

    fig = px.scatter(toPlot, x= x_axis, y = y_axis, trendline = 'ols',
                     trendline_color_override = linecolor, color_discrete_sequence = ['black'],
                     )
    fig.update_layout(hovermode=False)
    fig.update_layout(scene = go.layout.Scene(aspectratio = {'x':1, 'y':1}),
                      height = 400,
                      plot_bgcolor='rgba(0,0,0,0)')
    scatter_div = plot(fig, output_type="div",)

    context: dict = {'scatter_plot': scatter_div, 'title': title}

    return render(request, 'correlation_scatter.html', context)






class PathwayListView(generic.ListView):
    model = Pathway

    def get_context_data(self, **kwargs):
        context = super(PathwayListView, self).get_context_data(**kwargs)
        context['forTable'] = pathwayAndCelltype.objects.all().order_by('-averageScore')[:20]
        return context

class PathwayDetailView(generic.DetailView):

    def get(self, request, *args, **kwargs):
        
        hspcType = getHSPCcookie(request)

        SI_scores = specificity_index[hspcType.upper()]

        correction = getCorrectionCookie(request)

        pathway = get_object_or_404(Pathway, pk = kwargs['pk'])
        pactsS = pathwayAndCelltype.objects.filter(pathway= pathway, sorr = 's', hscPercent__gt=0.05, hspc_type = hspcType)
        pactsR = pathwayAndCelltype.objects.filter(pathway= pathway, sorr = 't', hscPercent__gt=0.05, hspc_type = hspcType)
        evidenceList = pathway.evidences
        pmids, keggs, pmcs = get_evidence_list(evidenceList)

        SI_sending = "Not calculated"
        SI_receiving = "Not calculated"
        if 'sending ' + pathway.name in SI_scores.index:
            SI_sending = SI_scores.loc['sending ' + pathway.name,]
        if 'receiving ' + pathway.name in SI_scores.index:
            SI_receiving = SI_scores.loc['receiving ' + pathway.name,]

        correlatedPathways1 = pathwayCorrelations.objects.filter(pathway1 = pathway, hspc_type = hspcType)
        correlatedPathways2 = pathwayCorrelations.objects.filter(pathway2 = pathway, hspc_type = hspcType)
        
        ligands = pathway.ligands.all()
        receptors = pathway.receptors.all()

        rDF = receiving_pathways_dictionary_db[pathway.name]
        sDF = sending_pathways_dictionary_db[pathway.name]

        fig_total = make_subplots(rows = 1, cols = 2,
                                  subplot_titles = ("HSPCs Receiving Signals (population)", "HSPCs Sending Signals (population)"),
                                  
                                  horizontal_spacing=0.5)

        fig_total.add_trace(go.Heatmap(z = rDF,
                                      x = rDF.columns,
                                      y = rDF.index,
                                      colorscale = 'Reds',
                                      zmid = 0.5,
                                      xgap = 1,
                                      ygap = 1,
                                      hovertemplate = "%{y}→%{x}<br>Normalized Interaction Potential: %{z}<extra></extra>",
                                      colorbar_x=0.25),
                                      row = 1, col = 1)
        #rheatmap_div = plot(rfig, output_type="div",)
        fig_total.add_trace(go.Heatmap(z = sDF,
                                      x = sDF.columns,
                                      y = sDF.index,
                                      colorscale = 'Blues',
                                      zmid = 0.5,
                                      xgap = 1,
                                      ygap = 1,
                                      hovertemplate = "%{y}←%{x}<br>Normalized Interaction Potential: %{z}<extra></extra>"),
                                      row = 1, col = 2)
        #sheatmap_div = plot(sfig, output_type="div",)
        fig_total = plot(fig_total, output_type="div",)
        if len(pactsS) > 0 or len(pactsR) > 0:
            plot_div = make_net_graph_JSON(pactsS, pactsR, hspcType, 'cts',correction, signaling_type = pathway.interaction_type)
            #plot_div = make_net_graph_spread(pactsS, pactsR, 'cts')
            context = {'pathway': pathway, 'pactsS': pactsS, 'pactsR': pactsR,
                        'plot_div': plot_div, 'ligands':ligands, 'receptors': receptors,
                        'pmids' : pmids, 'keggs' : keggs, 'pmcs' : pmcs,
                        'correlations1' : correlatedPathways1, 'correlations2' : correlatedPathways2, 'hspcType' : hspcType.upper(),
                        'rheatmap_div':fig_total, 'sheatmap_div':fig_total,
                        'SI_sending': SI_sending, 'SI_receiving': SI_receiving}
            return render(request, 'interactions/pathway_detail.html', context)
        
        else:
            return render(request, 'interactions/pathway_detail.html', context = {'pathway': pathway,'ligands':ligands, 'receptors': receptors, 'hspcType' : hspcType.upper(), 'SI_sending': SI_sending, 'SI_receiving': SI_receiving})

class CellClassListView(generic.ListView):
    model = cellClas
    
    def get_context_data(self, **kwargs):
        context = super(CellClassListView, self).get_context_data(**kwargs)
        context['forTable'] = pathwayAndCelltype.objects.all().order_by('-averageScore')[:20]
        return context

class CellClassDetailView(generic.DetailView):

    def get(self, request, *args, **kwargs):
        
        hspcType = getHSPCcookie(request)

        cellclass = get_object_or_404(cellClas, pk = kwargs['pk'])

        ctp = 'NA'
        if str(cellclass) in cell_type_proportions.keys():
            ctp = cell_type_proportions[str(cellclass)]
        phenoMeasurement = 'NA'
        inicheMeasurement = 'NA'
        if hspcType in ['hsc' ,'mpp', 'mep', 'gmp', 'cmp', 'clp']:
            if str(cellclass) in pheno_dicts[hspcType].keys():
                phenoMeasurement = pheno_dicts[hspcType][str(cellclass)]
                inicheMeasurement = iniche_dicts[hspcType][str(cellclass)]

        pactsS = pathwayAndCelltype.objects.filter(celltype= cellclass, sorr = 's', hscPercent__range=(0.05,2), hspc_type = hspcType)
        pactsR = pathwayAndCelltype.objects.filter(celltype= cellclass, sorr = 't', hscPercent__range=(0.05,2), hspc_type = hspcType)

        receiving_celltypes_dictionary_db

        celltype_name = cellclass.name
        if celltype_name == 'Osteoclasts':
            celltype_name = 'Osteo'

        rDF = receiving_celltypes_dictionary_db[celltype_name]
        sDF = sending_celltypes_dictionary_db[celltype_name]

        fig_total = make_subplots(rows = 1, cols = 2,
                                  subplot_titles = ("HSPCs Receiving Signals (population)", "HSPCs Sending Signals (population)"),
                                  
                                  horizontal_spacing=0.5)

        fig_total.add_trace(go.Heatmap(z = rDF,
                                      x = rDF.columns,
                                      y = rDF.index,
                                      colorscale = 'Reds',
                                      zmid = 0.5,
                                      xgap = 1,
                                      ygap = 1,
                                      hovertemplate = "%{y}→%{x}<br>Normalized Interaction Potential: %{z}<extra></extra>",
                                      colorbar_x=0.25),
                                      row = 1, col = 1)
        #rheatmap_div = plot(rfig, output_type="div",)
        fig_total.add_trace(go.Heatmap(z = sDF,
                                      x = sDF.columns,
                                      y = sDF.index,
                                      colorscale = 'Blues',
                                      zmid = 0.5,
                                      xgap = 1,
                                      ygap = 1,
                                      hovertemplate = "%{y}←%{x}<br>Normalized Interaction Potential: %{z}<extra></extra>"),
                                      row = 1, col = 2)
        #sheatmap_div = plot(sfig, output_type="div",)
        fig_total = plot(fig_total, output_type="div",)

        if len(pactsS) > 0 or len(pactsR) >0:
            plot_div = make_net_graph_JSON(pactsS, pactsR, hspcType, 'paths')
            #plot_div = make_net_graph_spread(pactsS, pactsR, 'paths')
            context = {'cellclass': cellclass, 'pactsS': pactsS, 'pactsR': pactsR,
                        'plot_div': plot_div, 'hspcType' : hspcType.upper(), 'cellTypeProp': ctp, 'phenoPvalue': phenoMeasurement, 'inicheValue' : inicheMeasurement,
                        'rheatmap_div':fig_total, 'sheatmap_div':fig_total}
            return render(request, 'interactions/cellclas_detail.html', context)
        
        else:
            return render(request, 'interactions/cellclas_detail.html', context = {'pathway': cellclass, 'hspcType' : hspcType.upper()})

class LigandDetailView(generic.DetailView):

    def get(self, request, *args, **kwargs):

        hspcType = getHSPCcookie(request)

        genename = get_object_or_404(Ligand, pk=kwargs['pk'])
        pathways = genename.pathways.all()
        paired = genename.receptors.all()
        evidenceList = genename.evidences
        pmids, keggs, pmcs = get_evidence_list(evidenceList)
        genemap, paths = make_gene_map_specific(pathways)

        gene_expression_link = "https://gexc.riken.jp/models/3/genes/" + str.lower(genename.name)

        context = {"genename": genename, 'pathways': pathways, 'paired': paired, 'link': gene_expression_link,
                    'lorr':'l', 'plot_div':genemap, 'paths':paths,
                    'pmids' : pmids, 'keggs' : keggs, 'pmcs' : pmcs}

        return render(request, 'interactions/gene_detail.html', context = context)

class ReceptorDetailView(generic.DetailView):

    def get(self, request, *args, **kwargs):

        genename = get_object_or_404(Receptor, pk=kwargs['pk'])
        pathways = genename.pathways.all()
        paired = genename.ligands.all()
        evidenceList = genename.evidences
        pmids, keggs, pmcs = get_evidence_list(evidenceList)
        genemap, paths = make_gene_map_specific(pathways)

        gene_expression_link = "https://gexc.riken.jp/models/3/genes/" + str.lower(genename.name)

        context = {"genename": genename, 'pathways': pathways, 'paired': paired, 'link': gene_expression_link,
                    'lorr':'r', 'plot_div':genemap, 'paths':paths,
                    'pmids' : pmids, 'keggs' : keggs, 'pmcs' : pmcs}

        return render(request, 'interactions/gene_detail.html', context = context)



MIN_SIZE = 20
MAX_SIZE = 50

def make_net_graph(sending, receiving, ct_or_p = 'cts'):
    """ 
    View demonstrating how to display a graph object
    on a web page with Plotly. 
    """
    g = nx.Graph()

    sending_dict_avgscore_hscP = {}
    receiving_dict_avgscore_hscP = {}

    hscPs = []

    edge_weights = []
    
    title_name = ""
    center_cell = ""

    for ct in sending:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = "HSCs"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type

        g.add_edge(center_cell,cn+"_S")
        sending_dict_avgscore_hscP[cn] = [ct.averageScore]
        sending_dict_avgscore_hscP[cn].append(ct.hscPercent)
        sending_dict_avgscore_hscP[cn].append(t)
        hscPs.append(ct.hscPercent)
        edge_weights.append(ct.averageScore)
        

    for ct in receiving:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = "HSCs"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type

        g.add_edge(cn+"_R", center_cell)
        receiving_dict_avgscore_hscP[cn] = [ct.averageScore]
        receiving_dict_avgscore_hscP[cn].append(ct.hscPercent)
        receiving_dict_avgscore_hscP[cn].append(t)
        hscPs.append(ct.hscPercent)
        edge_weights.append(ct.averageScore)

    minP = min(hscPs)
    maxP = max(hscPs)
    ewMin = min(edge_weights)
    ewMax = max(edge_weights) 
    #edge_weights = [(i-ewMin)*(9)/(ewMax-ewMin) + 1 for i in edge_weights]
    x, y = circle(len(g.nodes()))
    
    pos = {}
    k = 0
    for i in g.nodes():
        if i != center_cell:
            pos[i] = [x[k], y[k]]
            k += 1
        else:
            pos[i] = [0,0]
        

    #pos = nx.planar_layout(g)

    edge_x = []
    edge_y = []
    edge_colors = []
    check = "_R" if ct_or_p == "cts" else "_S"
    for edge in g.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]


        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

        if check in edge[1]:
            edge_colors.append("orchid")
        else:
            edge_colors.append("aqua")
        '''
        if len(g.edges()) > 2:
            if check in edge[1]:
                edge_x, edge_y = addEdge([x1,y1], [x0,y0], edge_x, edge_y, .8, 'end', .05, 30, 55)
                edge_colors.append("orchid")
            else:
                edge_x, edge_y = addEdge([x0,y0], [x1,y1], edge_x, edge_y, .8, 'end', .05, 30, 55)
                edge_colors.append("aqua")
        else:
            edge_x, edge_y = addEdge([x0,y0], [x1,y1], edge_x, edge_y, .95, None, .05, 30, 55)

    edge_trace = {}

    for i in range(len(edge_weights)):
        edge_trace['trace_' + str(i)]= go.Scatter(x = edge_x[i*9:i*9+9],
                                                  y = edge_y[i*9:i*9+9],
                                                  line=dict(width=2, color = edge_colors[i]),
                                                  hoverinfo='none',
                                                  mode='lines')
    '''
    edge_trace = {}
    for i in range(len(edge_weights)):
        edge_trace['trace_' + str(i)]= go.Scatter(x = edge_x[i*3:i*3+3],
                                                  y = edge_y[i*3:i*3+3],
                                                  line=dict(width=2, color = edge_colors[i]),
                                                  hoverinfo='none',
                                                  mode='lines')

    edge_trace = list(edge_trace.values())
    
    #edge_trace = go.Scatter(
    #x=edge_x, y=edge_y,
    #line=dict(width=2, color=edge_colors),
    #hoverinfo='none',
    #mode='lines')

    node_x, node_y = [], []
    for node in g.nodes():
        x,y = pos[node]
        node_x.append(x)
        node_y.append(y)

    marker_type = 'circle' if ct_or_p == "paths" else 'square'

    node_trace = go.Scatter(
        x = node_x, y = node_y,
        mode = 'markers',
        marker_symbol = marker_type,
        hoverinfo = 'text',
        marker=dict(
        showscale=True,
        colorscale='Greys',
        reversescale=False,
        color=[],
        line = dict(width = 5),
        size=20,
        colorbar=dict(
            thickness=15,
            title='Average Probability Score',
            xanchor='left',
            titleside='right'
        ),
    ))

    node_avgscore = []
    node_text = []
    node_sizes = []
    node_ecs = []
    for node in g.nodes:
        if node == center_cell:
            node_avgscore.append(max(edge_weights))
            node_sizes.append(50)
            node_text.append(center_cell)
            node_ecs.append('black')
        else:
            cell_name = node.split("_")[0]
            if "_S" in node:
                avgScore = sending_dict_avgscore_hscP[cell_name][0]
                avgScore_txt = "Average Score: {:.2f}".format(avgScore)
                hscP = sending_dict_avgscore_hscP[cell_name][1]
                hscP_txt = "HSC%: {:.2f}".format(hscP*100)
                node_avgscore.append(avgScore)
                node_sizes.append(normalize_dot_size(MIN_SIZE, MAX_SIZE, hscP, minP, maxP))
                node_text.append("<br>".join([cell_name, avgScore_txt, hscP_txt]))
                node_ecs.append(get_ec(sending_dict_avgscore_hscP[cell_name][2], ct_or_p))

            else:
                avgScore = receiving_dict_avgscore_hscP[cell_name][0]
                avgScore_txt = "Average Score: {:.2f}".format(avgScore)
                hscP = receiving_dict_avgscore_hscP[cell_name][1]
                hscP_txt = "HSC%: {:.2f}".format(hscP*100)
                node_avgscore.append(avgScore)
                node_sizes.append(normalize_dot_size(MIN_SIZE, MAX_SIZE, hscP, minP, maxP))
                node_text.append("<br>".join([cell_name, avgScore_txt, hscP_txt]))
                node_ecs.append(get_ec(receiving_dict_avgscore_hscP[cell_name][2], ct_or_p))
    
    node_trace.marker.color = node_avgscore
    node_trace.text = node_text
    node_trace.marker.size = node_sizes
    node_trace.marker.line.color = node_ecs

    #normalize_dot_size(min_size, max_size, val, min_p, max_p)

    edge_trace.append(node_trace)
    #data = [edge_trace, node_trace]
    layout={
                'title': '<b>' + title_name + ' Pathway (Hover for info)</b>',
                'titlefont_size':16,
                'showlegend':False,
                'hovermode':'closest',
                'margin':dict(b=5,l=5,r=5,t=40),
                'xaxis':dict(showgrid=False, zeroline=False, showticklabels=False),
                'yaxis':dict(showgrid=False, zeroline=False, showticklabels=False),
                'width':700
                
    }
    
    
    # Getting HTML needed to render the plot.
    plot_div = plot({'data': edge_trace, 'layout': layout}, 
                    output_type='div')

    return plot_div

def genes(request):
    ligand_list = Ligand.objects.all()
    receptor_list = Receptor.objects.all()
    genemap, paths = make_gene_map(ligand_list)
    context = {'ligand_list': ligand_list, 'receptor_list': receptor_list, 'plot_div':genemap, 'paths':paths}
    
    
    return render(request, 'genes.html', context= context)

def make_gene_map(ligand_list):
    
    ligands= []
    receptors = []
    pathway = []
    l_urls = []
    r_urls = []

    for l in ligand_list:
        for r in l.receptors.all():
            for p in l.pathways.all():
                ligands.append(l.name)
                receptors.append(r.name)
                pathway.append(p.name)
                l_urls.append(l.get_absolute_url())
                r_urls.append(r.get_absolute_url())

    toret = json.dumps([ligands, receptors, pathway, list(set(pathway)), l_urls, r_urls])
    upaths = list(set(pathway))
    upaths.sort()
    return toret, upaths

def make_gene_map_specific(pathway_list):
    
    ligands= []
    receptors = []
    pathway = []
    l_urls = []
    r_urls = []

    for p in pathway_list:
        for l in p.ligands.all():
            for r in l.receptors.all():
                ligands.append(l.name)
                receptors.append(r.name)
                pathway.append(p.name)
                l_urls.append(l.get_absolute_url())
                r_urls.append(r.get_absolute_url())

    toret = json.dumps([ligands, receptors, pathway, list(set(pathway)), l_urls, r_urls])
    upaths = list(set(pathway))
    upaths.sort()
    return toret, upaths

def make_net_graph_spread(sending, receiving, ct_or_p = 'cts'):
    """ 
    View demonstrating how to display a graph object
    on a web page with Plotly. 
    """
    g = nx.Graph()

    sending_dict_avgscore_hscP = {}
    receiving_dict_avgscore_hscP = {}

    hscPs = []

    edge_weights = []
    
    title_name = ""
    center_cell = ""

    for ct in sending:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = "HSCs"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type
        g.add_node(center_cell)
        g.add_node(cn+"_S")
        sending_dict_avgscore_hscP[cn] = [ct.averageScore]
        sending_dict_avgscore_hscP[cn].append(ct.hscPercent)
        sending_dict_avgscore_hscP[cn].append(t)
        hscPs.append(ct.hscPercent)
        edge_weights.append(ct.averageScore)
        
    for ct in receiving:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = "HSCs"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type
        g.add_node(center_cell)
        g.add_node(cn+"_R")
        receiving_dict_avgscore_hscP[cn] = [ct.averageScore]
        receiving_dict_avgscore_hscP[cn].append(ct.hscPercent)
        receiving_dict_avgscore_hscP[cn].append(t)
        hscPs.append(ct.hscPercent)
        edge_weights.append(ct.averageScore)
    
    minP = min(hscPs)
    maxP = max(hscPs)

    x, y = circle(len(g.nodes())-1)
    
    pos = {}
    k = 0
    for i in g.nodes():
        if i != center_cell:
            pos[i] = [x[k], y[k]]
            k += 1
        else:
            pos[i] = [0,0]


    for ct in sending:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = "HSCs"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type

        g.add_edge(center_cell,cn+"_S")        

    for ct in receiving:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = "HSCs"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type

        g.add_edge(cn+"_R", center_cell)


    edge_x = []
    edge_y = []
    edge_colors = []
    check = "_R" if ct_or_p == "cts" else "_S"
    for edge in g.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]


        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

        if check in edge[1]:
            edge_colors.append("orchid")
        else:
            edge_colors.append("aqua")
        '''
        if len(g.edges()) > 2:
            if check in edge[1]:
                edge_x, edge_y = addEdge([x1,y1], [x0,y0], edge_x, edge_y, .8, 'end', .05, 30, 55)
                edge_colors.append("orchid")
            else:
                edge_x, edge_y = addEdge([x0,y0], [x1,y1], edge_x, edge_y, .8, 'end', .05, 30, 55)
                edge_colors.append("aqua")
        else:
            edge_x, edge_y = addEdge([x0,y0], [x1,y1], edge_x, edge_y, .95, None, .05, 30, 55)

    edge_trace = {}

    for i in range(len(edge_weights)):
        edge_trace['trace_' + str(i)]= go.Scatter(x = edge_x[i*9:i*9+9],
                                                  y = edge_y[i*9:i*9+9],
                                                  line=dict(width=2, color = edge_colors[i]),
                                                  hoverinfo='none',
                                                  mode='lines')
    '''
    edge_trace = {}
    for i in range(len(edge_weights)):
        edge_trace['trace_' + str(i)]= go.Scatter(x = edge_x[i*3:i*3+3],
                                                  y = edge_y[i*3:i*3+3],
                                                  line=dict(width=2, color = edge_colors[i]),
                                                  hoverinfo='none',
                                                  mode='lines')

    edge_trace = list(edge_trace.values())
    
    #edge_trace = go.Scatter(
    #x=edge_x, y=edge_y,
    #line=dict(width=2, color=edge_colors),
    #hoverinfo='none',
    #mode='lines')

    node_x, node_y = [], []
    for node in g.nodes():
        x,y = pos[node]
        node_x.append(x)
        node_y.append(y)

    marker_type = 'circle' if ct_or_p == "paths" else 'diamond'
    marker_labels = [i[:-2] for i in list(g.nodes())[1:]]
    marker_labels.insert(0, "")

    node_trace = go.Scatter(
        x = node_x, y = node_y,
        mode = 'markers+text',
        marker_symbol = marker_type,
        text = marker_labels,
        textposition="bottom center",
        hoverinfo = 'text',
        marker=dict(
        showscale=True,
        colorscale='Greys',
        reversescale=False,
        color=[],
        line = dict(width = 5),
        size=20,
        opacity = 1,
        colorbar=dict(
            thickness=15,
            title='Average Probability Score',
            xanchor='left',
            titleside='right'
        ),
    ))

    node_avgscore = []
    node_text = []
    node_sizes = []
    node_ecs = []
    for node in g.nodes:
        if node == center_cell:
            node_avgscore.append(max(edge_weights))
            node_sizes.append(50)
            node_text.append(center_cell)
            node_ecs.append('black')
        else:
            cell_name = node.split("_")[0]
            if "_S" in node:
                avgScore = sending_dict_avgscore_hscP[cell_name][0]
                avgScore_txt = "Average Score: {:.2f}".format(avgScore)
                hscP = sending_dict_avgscore_hscP[cell_name][1]
                hscP_txt = "HSC%: {:.2f}".format(hscP*100)
                node_avgscore.append(avgScore)
                node_sizes.append(normalize_dot_size(MIN_SIZE, MAX_SIZE, hscP, minP, maxP))
                node_text.append("<br>".join([cell_name, avgScore_txt, hscP_txt]))
                node_ecs.append(get_ec(sending_dict_avgscore_hscP[cell_name][2], ct_or_p))

            else:
                avgScore = receiving_dict_avgscore_hscP[cell_name][0]
                avgScore_txt = "Average Score: {:.2f}".format(avgScore)
                hscP = receiving_dict_avgscore_hscP[cell_name][1]
                hscP_txt = "HSC%: {:.2f}".format(hscP*100)
                node_avgscore.append(avgScore)
                node_sizes.append(normalize_dot_size(MIN_SIZE, MAX_SIZE, hscP, minP, maxP))
                node_text.append("<br>".join([cell_name, avgScore_txt, hscP_txt]))
                node_ecs.append(get_ec(receiving_dict_avgscore_hscP[cell_name][2], ct_or_p))
    
    node_trace.marker.color = node_avgscore
    node_trace.hovertext = node_text
    node_trace.marker.size = node_sizes
    node_trace.marker.line.color = node_ecs

    #normalize_dot_size(min_size, max_size, val, min_p, max_p)

    edge_trace.append(node_trace)
    #data = [edge_trace, node_trace]
    layout={
                'title': '<b>' + title_name + ' Pathway (Hover for info)</b>',
                'titlefont_size':16,
                'showlegend':False,
                'hovermode':'closest',
                'margin':dict(b=5,l=5,r=5,t=40),
                'xaxis':dict(showgrid=False, zeroline=False, showticklabels=False),
                'yaxis':dict(showgrid=False, zeroline=False, showticklabels=False),
                'width':700
                
    }
    
    
    # Getting HTML needed to render the plot.
    plot_div = plot({'data': edge_trace, 'layout': layout}, 
                    output_type='div')

    return plot_div

def make_net_graph_corr(correlations):

    p1names = []
    p1colors = []
    p1ids = []
    p1urls = []

    p2names = []
    p2colors = []
    p2ids = []
    p2urls = []

    corrs = []
    pvals = []
    edge_colors = []    

    p1names_n = []
    p1colors_n = []
    p1ids_n = []
    p1urls_n = []

    p2names_n = []
    p2colors_n = []
    p2ids_n = []
    p2urls_n = []

    corrs_n = []
    pvals_n = []
    edge_colors_n = []    

    for corr in correlations:
        pathway1_name = corr.pathway1.name
        pathway1_id = pathway1_name + "_" + corr.p1a.upper()
        pathway1_color = "#9e4d9e" if corr.p1a.upper() == "S" else "#458045"
        pathway2_name = corr.pathway2.name
        pathway2_id = pathway2_name + "_" + corr.p2a.upper()
        pathway2_color = "#9e4d9e" if corr.p2a.upper() == "S" else "#458045"
        corr_val = corr.correlation
        p_val = corr.pval
        edge_color = "red" if corr_val > 0 else "blue"

        if corr_val > 0:
            p1names.append(pathway1_name)
            p1colors.append(pathway1_color)
            p1ids.append(pathway1_id)
            p1urls.append(corr.pathway1.get_absolute_url)

            p2names.append(pathway2_name)
            p2colors.append(pathway2_color)
            p2ids.append(pathway2_id)
            p2urls.append(corr.pathway2.get_absolute_url)

            corrs.append(corr_val)
            pvals.append(p_val)
            edge_colors.append(edge_color)
        
        else:
            p1names_n.append(pathway1_name)
            p1colors_n.append(pathway1_color)
            p1ids_n.append(pathway1_id)
            p1urls_n.append(corr.pathway1.get_absolute_url)

            p2names_n.append(pathway2_name)
            p2colors_n.append(pathway2_color)
            p2ids_n.append(pathway2_id)
            p2urls_n.append(corr.pathway2.get_absolute_url)

            corrs_n.append(corr_val)
            pvals_n.append(p_val)
            edge_colors_n.append(edge_color)
    
    corrs, p1names, p1colors, p1ids, p1urls, p2names, p2colors, p2ids, p2urls,  pvals, edge_colors = zip(*sorted(zip(corrs, p1names, p1colors, p1ids, p1urls, p2names, p2colors, p2ids, p2urls, pvals, edge_colors), reverse=True))

    corrs_n, p1names_n, p1colors_n, p1ids_n, p1urls_n, p2names_n, p2colors_n, p2ids_n, p2urls_n, pvals_n, edge_colors_n = zip(*sorted(zip(corrs_n, p1names_n, p1colors_n, p1ids_n, p1urls_n, p2names_n, p2colors_n, p2ids_n, p2urls_n, pvals_n, edge_colors_n)))

    toret = json.dumps([list(x) for x in [corrs, p1names, p1colors, p1ids, p2names, p2colors, p2ids, pvals, edge_colors, corrs_n, p1names_n, p1colors_n, p1ids_n, p2names_n, p2colors_n, p2ids_n, pvals_n, edge_colors_n, p1urls, p2urls, p1urls_n, p2urls_n]])

    return toret

def make_net_graph_JSON(sending, receiving, hspc_type, ct_or_p = 'cts', correction = 'none', signaling_type = 'c'):

    max_score = 0
    min_score = 100

    max_score_facs = 0
    min_score_facs = 100

    max_score_pheno = 0
    min_score_pheno = 0

    node_names = []
    to_names = []
    from_names = []

    sending_dict_avgscore_hscP = {}
    receiving_dict_avgscore_hscP = {}

    hscPs = []

    edge_weights = []
    
    title_name = ""
    center_cell = ""

    for ct in sending:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = hspc_type.upper() + "s"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type

        if ct.averageScore > max_score:
            max_score = ct.averageScore
        if ct.averageScore < min_score:
            min_score = ct.averageScore

        if ct_or_p == 'cts':
            if cn in cell_type_proportions.keys():
                if ct.averageScore*cell_type_proportions[cn] > max_score_facs:
                    max_score_facs = ct.averageScore*cell_type_proportions[cn]
                if ct.averageScore*cell_type_proportions[cn] < min_score_facs:
                    min_score_facs = ct.averageScore*cell_type_proportions[cn]
            if hspc_type in ['hsc', 'mpp', 'mep', 'gmp', 'cmp', 'clp']:
                if signaling_type == 'c':
                    pheno_dict_to_use = iniche_dicts[hspc_type]
                if signaling_type == 's':
                    pheno_dict_to_use = pheno_dicts[hspc_type]
                if signaling_type == 'e':
                    pheno_dict_to_use = pheno_dicts[hspc_type]
                if cn in pheno_dict_to_use.keys() and cn in cell_type_proportions.keys():
                    if ct.averageScore*cell_type_proportions[cn]*(pheno_dict_to_use[cn]) > max_score_pheno:
                        max_score_pheno = ct.averageScore*cell_type_proportions[cn]*(pheno_dict_to_use[cn])
                    if ct.averageScore*cell_type_proportions[cn]*(pheno_dict_to_use[cn]) < min_score_pheno:
                        min_score_pheno = ct.averageScore*cell_type_proportions[cn]*(pheno_dict_to_use[cn])

        sending_dict_avgscore_hscP[cn] = [ct.averageScore]
        sending_dict_avgscore_hscP[cn].append(ct.hscPercent)
        sending_dict_avgscore_hscP[cn].append(t)
        if ct_or_p == 'cts':
            sending_dict_avgscore_hscP[cn].append(ct.celltype.get_absolute_url)
        else:
            sending_dict_avgscore_hscP[cn].append(ct.pathway.get_absolute_url)

        hscPs.append(ct.hscPercent)
        edge_weights.append(ct.averageScore)
        
    for ct in receiving:
        if ct_or_p == 'cts':
            cn = ct.celltype.name
            title_name = ct.pathway.name
            center_cell = hspc_type.upper() + "s"
            t = ct.celltype.cell_type
        else:
            cn = ct.pathway.name
            title_name = ct.celltype.name
            center_cell = title_name
            t = ct.pathway.interaction_type

        if ct.averageScore > max_score:
            max_score = ct.averageScore
        if ct.averageScore < min_score:
            min_score = ct.averageScore

        if ct_or_p == 'cts':
            if cn in cell_type_proportions.keys():
                if ct.averageScore*cell_type_proportions[cn] > max_score_facs:
                    max_score_facs = ct.averageScore*cell_type_proportions[cn]
                if ct.averageScore*cell_type_proportions[cn] < min_score_facs:
                    min_score_facs = ct.averageScore*cell_type_proportions[cn]
            if hspc_type in ['hsc', 'mpp', 'mep', 'gmp', 'cmp', 'clp']:
                if signaling_type == 'c':
                    pheno_dict_to_use = iniche_dicts[hspc_type]
                if signaling_type == 's':
                    pheno_dict_to_use = pheno_dicts[hspc_type]
                if signaling_type == 'e':
                    pheno_dict_to_use = pheno_dicts[hspc_type]
                if cn in pheno_dict_to_use.keys() and cn in cell_type_proportions.keys():
                    if ct.averageScore*cell_type_proportions[cn]*(pheno_dict_to_use[cn]) > max_score_pheno:
                        max_score_pheno = ct.averageScore*cell_type_proportions[cn]*(pheno_dict_to_use[cn])
                    if ct.averageScore*cell_type_proportions[cn]*(pheno_dict_to_use[cn]) < min_score_pheno:
                        min_score_pheno = ct.averageScore*cell_type_proportions[cn]*(pheno_dict_to_use[cn])

        receiving_dict_avgscore_hscP[cn] = [ct.averageScore]
        receiving_dict_avgscore_hscP[cn].append(ct.hscPercent)
        receiving_dict_avgscore_hscP[cn].append(t)
        if ct_or_p == 'cts':
            receiving_dict_avgscore_hscP[cn].append(ct.celltype.get_absolute_url)
        else:
            receiving_dict_avgscore_hscP[cn].append(ct.pathway.get_absolute_url)
        
        hscPs.append(ct.hscPercent)
        edge_weights.append(ct.averageScore)


    node_names = list(set(node_names))
    
    node_names = []
    fill_color = []
    avgscore = []
    edges = []
    hscPercent = []
    height = []
    ecs = []
    urlList = []
    k = 0


    if signaling_type == 'c':
        pheno_dict_to_use = iniche_dicts[hspc_type]
    if signaling_type == 's':
        pheno_dict_to_use = pheno_dicts[hspc_type]
    if signaling_type == 'e':
        pheno_dict_to_use = pheno_dicts[hspc_type]

    for i in sending_dict_avgscore_hscP.keys():

        correction_coeff = 1
        
        if correction == 'phenocycler':
            if hspc_type not in pheno_dicts.keys() or i not in cell_type_proportions:
                break
            if i in KLS_spatial_pvalue.keys():
                correction_coeff = cell_type_proportions[i]*(pheno_dict_to_use[i])

                range = max_score_pheno - min_score_pheno
                a = sending_dict_avgscore_hscP[i][0]*correction_coeff*(pheno_dict_to_use[i])
                if range != 0:
                    a = (a - min_score) / range
                range2 = 9999
                a = (a * range2)
                fill_color.append(mpl.colors.to_hex(greys(int(a))))
            else:
                continue
        else:
            range = max_score - min_score
            a = sending_dict_avgscore_hscP[i][0]
            if range != 0:
                a = (a - min_score) / range
            range2 = 9999
            a = (a * range2)
            fill_color.append(mpl.colors.to_hex(greys(int(a))))


        node_names.append(i + "_" + str(k))
        ecs.append(get_ec(sending_dict_avgscore_hscP[i][2], ct_or_p))
        height.append(sending_dict_avgscore_hscP[i][1]*60 + 10)
        
        avgscore.append(sending_dict_avgscore_hscP[i][0]*correction_coeff)
        hscPercent.append(100*sending_dict_avgscore_hscP[i][1])
        edges.append('#8FAADC')
        to_names.append(i +"_" + str(k))
        from_names.append(center_cell)
        urlList.append(sending_dict_avgscore_hscP[i][3])
        k += 1

    for i in receiving_dict_avgscore_hscP.keys():

        correction_coeff = 1

        if correction == 'phenocycler':
            if hspc_type not in pheno_dicts.keys() or i not in cell_type_proportions:
                break
            if i in KLS_spatial_pvalue.keys() and i in cell_type_proportions.keys():
                correction_coeff = cell_type_proportions[i]*(pheno_dict_to_use[i])

                range = max_score_pheno - min_score_pheno
                a = receiving_dict_avgscore_hscP[i][0]*correction_coeff*(pheno_dict_to_use[i])
                if range != 0:
                    a = (a - min_score) / range
                range2 = 9999
                a = (a * range2)
                fill_color.append(mpl.colors.to_hex(greys(int(a))))
            else:
                continue
        else:
            range = max_score - min_score
            a = receiving_dict_avgscore_hscP[i][0]
            if range != 0:
                a = (a - min_score) / range
            range2 = 9999
            a = (a * range2)
            fill_color.append(mpl.colors.to_hex(greys(int(a))))


        node_names.append(i +"_" + str(k))
        ecs.append(get_ec(receiving_dict_avgscore_hscP[i][2], ct_or_p))
        height.append(receiving_dict_avgscore_hscP[i][1]*60 + 10)

        avgscore.append(receiving_dict_avgscore_hscP[i][0]*correction_coeff)
        hscPercent.append(100*receiving_dict_avgscore_hscP[i][1])
        from_names.append(i +"_" + str(k))
        edges.append('#FF8181')
        to_names.append(center_cell)
        urlList.append(receiving_dict_avgscore_hscP[i][3])
        k += 1
    node_names.append(center_cell)
    ecs.append('3 black')
    fill_color.append("black")
    avgscore.append("NaN")
    hscPercent.append("NaN")
    height.append(70)
    urlList.append("")

    toret = json.dumps([node_names, to_names, from_names, ecs, fill_color, avgscore, hscPercent, edges, height, urlList])

    return toret