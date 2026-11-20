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

def plotSpc1DBurgers(iPlot:int, dfData:pd.DataFrame, batchName:str, varName:str):
    # Control
    data = []; N = 360;
    vPi = math.pi; xVec = np.arange(0, (1 + 1/N) * 2 * vPi, 2 * vPi / N)

    # Data
    plotCols= sorted([col for col in dfData.columns if str(col).strip().isdigit()], key=lambda c: int(str(c).strip()))

    # Calculate
    for i, row in dfData.iterrows():
        # Data
        mData = row.loc[plotCols].dropna()
        uHat = mData.astype(str).astype(complex).to_numpy()
        Nk = len(uHat) - 1

        # Calculate
        uData = np.full_like(xVec, uHat[0].real, dtype=np.float64)
        for k in range(1, Nk):
            kx = k * xVec
            uData += 2 * (uHat[k].real * np.cos(kx) - uHat[k].imag * np.sin(kx))

        # Cleanup
        data.append({varName: row[varName], 'Solver': row['Solver'], 'label': row['label'], 'Time': row['Time'], 'uData': uData})

    # DataFrame
    dfPlot = pd.DataFrame(data)

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
    for i, row in dfPlot.iterrows():
        # Data
        if varName != "N":
            if row['Solver'] == "DNS": NFixed = 100
            if row['Solver'] == "LES": NFixed = 20
        vVal = int(row[varName]) if varName == 'N' else int(NFixed);

        # Plot
        style = stylePlt.get(row['Solver'], stylePlt['DEFAULT']); stylePlt[row['Solver']]['bData'] = True
        if float(row['Time']) > 9.99:
            tempax = style['twin']; templinestyle = ':'
            templabel = f"{row['label']} (Diverged)"
            stylePlt[row['Solver']]['bTwin'] = True; 
        else:
            tempax = style['plot']; templinestyle = style['linestyle']
            templabel = row['label']
        tempax.plot(xVec, row['uData'], label=templabel, marker=style['marker'], linestyle=templinestyle, ms=3, linewidth=1.5)

    # Cleanup
    for key, val in stylePlt.items():
        # Control
        if key == 'DEFAULT': continue
        handTwin, labTwin = [], []

        # Configuration (Twin)
        if val['bTwin']:
            val['twin'].set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
            val['twin'].set_xticklabels([r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
            val['twin'].tick_params(axis='both', labelsize=PG.iFontSizeTicks)
            val['twin'].set_ylabel(f"$u(x)$ (Diverged)", fontsize=PG.iFontSizeLabel)
            handTwin, labTwin = val['twin'].get_legend_handles_labels()
        else:
            val['twin'].yaxis.set_visible(False)
        
        # Plot Toggle
        if not val['bData']: 
            val['plot'].set_visible(False)
            val['twin'].set_visible(False)

        # Configuration (Plot)
        val['plot'].set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi]); 
        val['plot'].set_xticklabels([r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
        val['plot'].tick_params(axis='both', labelsize=PG.iFontSizeTicks)
        val['plot'].set_xlabel('$x$', fontsize=PG.iFontSizeLabel); val['plot'].set_ylabel(f"$u(x)$", fontsize=PG.iFontSizeLabel)
        val['plot'].grid(True, which="both", ls="--", alpha=0.5)
        plt.tight_layout()

        # Legend
        handPlot, labPlot = val['plot'].get_legend_handles_labels()
        val['plot'].legend(handPlot + handTwin, labPlot + labTwin, fontsize=PG.iFontSizeLegend);

    # Save Plot
    writePath = PG.pathPlt / "Burgers"; writePath.mkdir(parents=True, exist_ok=True)
    tempName = "Plot_" + batchName + "_" + str(iPlot) + ".png"; 
    if not os.path.exists(writePath / tempName):
        print("Exporting image ..."); plt.savefig(writePath / tempName); print(f"File saved to: {writePath / tempName}")


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
