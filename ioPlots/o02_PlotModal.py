### Imports
import os
import re
import sys
import math
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

### Self-Imports
import o09_PlotGlobal as PG

def plotRef1DBurgers(data:pd.DataFrame):
    # Receives DataFrame with data
    # Generates plots and stores them in directory
    quit()

def plotSpcPreProcess(dfData:pd.DataFrame, oData:list):

    pass

def plotSpc1DBurgers(iPlot:int, dfData:pd.DataFrame, batchName:str, varName:str):
    # Control
    data = [];
    vPi = math.pi; xVec = np.arange(0, (1 + 1/36) * 2 * vPi, 2 * vPi / 36)

    # LOOP OVER RANGE AND CREATE SPATIAL VELOCITY
    # THEN JUST PLOT THAT

    # This will plot the spatial velocity instead of the energy spectrum
    # xVec goes from [0, 2 * pi] -- Periodic BC
    # yVec is calculated from sum

    print('Spatial analysis.')

    for xVal in xVec:
        yTemp = [];
        # PENDING HERE, RECONSTRUCT SPATIAL VELOCITY

    # # Calculate
    # for i, row in dfData.iterrows():
        


def plotNum1DBurgers(iPlot:int, dfData:pd.DataFrame, batchName:str, varName:str):
    # Styles
    fig, (axDNS, axLES) = plt.subplots(1, 2, figsize=(18, 6), sharey=True)
    axDNS_t = axDNS.twinx(); axLES_t = axLES.twinx();
    stylePlt = {
        'DNS': {'plot': axDNS, 'twin': axDNS_t, 'bTwin': False, 'bData': False, 'marker': 'x', 'linestyle': '--'},
        'LES': {'plot': axLES, 'twin': axLES_t, 'bTwin': False, 'bData': False, 'marker': 'x', 'linestyle': '--'},
        'DEFAULT': {'marker': 'x', 'linestyle': '-'},
    }

    # Plotting
    NFixed = 0
    for i, row in dfData.iterrows():
        # Data
        if varName != "N":
            if row['Solver'] == "DNS": NFixed = 100
            if row['Solver'] == "LES": NFixed = 20
        vVal = int(row[varName]) if varName == 'N' else int(NFixed);
        xVec = np.arange(1, vVal+1); 
        yCol = [str(k) for k in range(1, vVal+1)]; yVec = row[yCol].values.astype(float)

        # Plot
        style = stylePlt.get(row['Solver'], stylePlt['DEFAULT']); stylePlt[row['Solver']]['bData'] = True
        if np.max(yVec) > 1:
            tempax = style['twin']; templinestyle = ':'
            templabel = f"{row['label']} (Diverged)"
            stylePlt[row['Solver']]['bTwin'] = True; 
        else:
            tempax = style['plot']; templinestyle = style['linestyle']
            templabel = row['label']
        tempax.plot(xVec, yVec, label=templabel, marker=style['marker'], linestyle=templinestyle, ms=3, linewidth=1.5)

    # Cleanup
    for key, val in stylePlt.items():
        # Control
        if key == 'DEFAULT': continue
        handTwin, labTwin = [], []

        # Configuration (Twin)
        if val['bTwin']:
            val['twin'].set_xscale('log'); val['twin'].set_yscale('log');
            val['twin'].tick_params(axis='both', labelsize=PG.iFontSizeTicks)
            val['twin'].set_ylabel(f"$E_k$ (Diverged)", fontsize=PG.iFontSizeLabel)
            handTwin, labTwin = val['twin'].get_legend_handles_labels()
        else:
            val['twin'].yaxis.set_visible(False)
        
        # Plot Toggle
        if not val['bData']: 
            val['plot'].set_visible(False)
            val['twin'].set_visible(False)

        # Slope
        varVal = 100
        xVec = np.arange(1, varVal+1); yVec = [1/int(k)/int(k) for k in xVec]
        val['plot'].plot(xVec, yVec, label='Slope m=2', linestyle='-.', linewidth=1.5)

        # Configuration (Plot)
        val['plot'].set_xscale('log'); val['plot'].set_yscale('log');
        val['plot'].tick_params(axis='both', labelsize=PG.iFontSizeTicks)
        val['plot'].set_xlabel('$k$', fontsize=PG.iFontSizeLabel); val['plot'].set_ylabel(f"$E_k$", fontsize=PG.iFontSizeLabel)
        val['plot'].grid(True, which="both", ls="--", alpha=0.5)
        plt.tight_layout()

        # Legend
        handPlot, labPlot = val['plot'].get_legend_handles_labels()
        val['plot'].legend(handPlot + handTwin, labPlot + labTwin, fontsize=PG.iFontSizeLegend);

    # Save Plot
    writePath = PG.pathPlt / "Burgers"; writePath.mkdir(parents=True, exist_ok=True)
    tempName = "Plot_" + batchName + ".png"; 
    if not os.path.exists(writePath / tempName):
        print("Exporting image ..."); plt.savefig(writePath / tempName); print(f"File saved to: {writePath / tempName}")
