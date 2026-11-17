### Imports
import os
import re
import sys
import numpy as np
import pandas as pd
import pathlib as pl

def importFullData(oPath:list, oData:list): # Import entire .csv
    # Import
    df = pd.read_csv(oPath)
    if df.empty: return 
    oData = df.to_dict()

    return oData 

def importLastRow(oDir:pl.PosixPath, oPath:pl.PosixPath, oData:list, varName:str): # Append last row from .csv
    # Import
    df = pd.read_csv(oPath)
    if df.empty: return
    lRow = df.iloc[-1].to_dict()
    
    # Data
    argVar = re.search(rf'{varName}(-?\d+(?:\.\d+)?)', oDir.name); lRow[varName] = float(argVar.group(1))
    argSolver = re.search(r'(DNS|LES)', oDir.name); lRow['Solver'] = argSolver.group(1)
    lRow['label'] = f"{lRow['Solver']} ({varName}={lRow[varName]})"
    oData.append(lRow)

    return oData

########## Main Plotter File ########## [OLD]
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

