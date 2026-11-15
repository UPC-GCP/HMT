### Imports
import os
import re
import sys
import numpy as np
import pandas as pd
from pathlib import Path

def importFullData(oPath:list, oData:list): # oList -- Import entire .csv
    # Import
    df = pd.read_csv(oPath)
    if df.empty: return 

    # Receives oData empty object
    # Returns it as full .csv file

    return oData 

def importLastRow(oPath:list, oData:list): # oList -- Append last row from .csv
    # Import
    df = pd.read_csv(oPath)
    if df.empty: return

    # Receives oData empty/incomplete object
    # Returns it with appended last row

    return oData

########## Main Plotter File ##########
### Use this file to generate the desired plots.
### Plot functions should be added in PlotLib.py

########## Types of Plots ##########
### Animation (Animate Map): Reads coordinates from file header
### Snapshot (Plot Map): If timeStep not specified, assumes last frame
### Profile (Plot Temperature Profile): If timeStep not specified assumes last frame, needs separate file with a single row of temperatures
### Point (Plot Point Evolution): Plot data from Plot_0_Point.

########## File Inputs ##########
### Path: Folder where files are located
### bPoint: Plot point data
### bAnimate: Create animation of temperature map evolution (Probe_1_Map)
### vSnapshot: Vector of instants for temperature map snapshot
### bVelocity: Create velocity field map and snapshot
### bNumerical: Plot numerical study
### sVar: Variable for numerical study

# Safety Check
# if len(sys.argv) != 5 and len(sys.argv) != 6 and len(sys.argv) != 8: print("Arguments passed incorrectly. (Expected values: Path, bPoint, bAnimate, vSnapshot)"); print("Additional arguments available. (Optional values: bVelocity, iNumerical, sVar)"); quit()
# if len(sys.argv) != 0: print("Arguments passed incorrectly.")

