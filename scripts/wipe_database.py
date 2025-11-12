from tkinter import E
from interactions.models import Pathway, Receptor, Ligand, cellClas, pathwayAndCelltype, pathwayCorrelations
import pandas as pd
import numpy as np
# import all the code to properly extract relevant data from cellchat outputs

def run():

    # Set up pathways and cell types and their pairs

    Pathway.objects.all().delete()
    cellClas.objects.all().delete()
    pathwayAndCelltype.objects.all().delete()
    Ligand.objects.all().delete()
    Receptor.objects.all().delete()
    pathwayCorrelations.objects.all().delete()
