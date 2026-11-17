### Imports
import os
import re
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

### Self-Imports
import o09_PlotGlobal as PG
import o01_PlotImport as PI
import o02_PlotModal as PM

##### Plotting #####
def newPlotReference(plotName:str): # 1D Burgers Reference Plots
    # Control
    data = []; fNames = ['*cfgBurgersMeshDNSN20', '*cfgBurgersMeshDNSN100', '*cfgBurgersKolCLESCK0.05', '*cfgBurgersKolCLESCK0.4223']
    pDir = [d for fName in fNames for d in PG.pathRes.glob(fName) if d.is_dir()]
    
    print(PG.pathRes, PG.pathPlt)

def newPlotNumerical(plotName:str, batchName:str, varName:str): # 1D Burgers Parametric Studies
    # Control
    data = []; iPlot = None;

    # Directories
    pDir = [d for d in PG.pathRes.glob(f"*{batchName}*") if d.is_dir()]
    for vDir in pDir:
        # Files
        pFiles = list(vDir.glob("*.csv"))
        for vFile in pFiles:
            # Control
            if vFile.name.find(plotName) == -1: continue
            iPlot = re.search(rf'Probe_(\d+)', vFile.name).group(1)

            # Data
            data = PI.importLastRow(vDir, vFile, data, varName)

    # DataFrame
    dfData = pd.DataFrame(data); dfData = dfData.sort_values(by=['Solver', varName]).reset_index(drop=True)

    # Plotting
    if iPlot == '0': PM.plotNum1DBurgers(iPlot, dfData, batchName, varName)
    elif iPlot == '1': PM.plotSpc1DBurgers(iPlot, dfData, batchName, varName)

##### Instructions #####
### Plotting Toolkit Structure
### If any modifications need to be made they can be included here and separated into specific modules
# o09_PlotGlobal: Global constants
# o01_PlotImport: importFullData(), importLastRow() -- Pass empty/incomplete list, returns with data
# o02_PlotModal: plotRef1DBurgers(), plotNum1DBurgers() -- Pass DataFrame, generates plots
# o02_PlotGeometry:
# o02_PlotAnimation:

### Input Arguments
# newPlotReference: plotName:str (e.g. Probe_0)
# newPlotNumerical: batchName:str (e.g. cfgBurgersMesh), varName:str (e.g. N)

##### Run #####
iArgs = len(sys.argv)
if iArgs == 1: print("Program Shutdown: Missing arguments.")
elif iArgs == 2: newPlotReference(sys.argv[1])
elif iArgs == 4: newPlotNumerical(sys.argv[1], sys.argv[2], sys.argv[3])

########## OLD CODE ##########
##### Use this as reference

def plotFinalBurgers(): # 1D Burgers Reference Plots
    # Control
    data = []
    # plotImport.import*(path, List) -- List object must be initialized here
    # Objects must be initialized a directed from here, should pass file directly instead of folder
    # Should know which option of files I have first and then loop over the alternatives

    # Import Data
    fNames = ['*cfgBurgersMeshDNSN20', '*cfgBurgersMeshDNSN100', '*cfgBurgersKolCLESCK0.05', '*cfgBurgersKolCLESCK0.4223']
    pDir = [d for fName in fNames for d in pathRes.glob(fName) if d.is_dir()]
    for vDir in pDir:
        # Find file
        pFiles = list(vDir.glob('*.csv'))
        for pFile in pFiles:
            # Import File
            if pFile.name.find("0") == -1: continue
            df = pd.read_csv(pFile)
            if df.empty: continue
            lRow = df.iloc[-1].to_dict()

            # Data
            argSolver = re.search(r'(DNS|LES)', vDir.name); lRow['Solver'] = argSolver.group(1)
            if lRow['Solver'] == "DNS": argN = re.search(rf"{lRow['Solver']}N(\d+)", vDir.name);
            lRow['N'] = float(argN.group(1) if lRow['Solver'] == "DNS" else 20)
            if lRow['Solver'] == "LES": argCk = re.search(rf"{lRow['Solver']}CK(\d+(?:\.\d+)?)", vDir.name); lRow['Ck'] = argCk.group(1)
            sLab = f"{lRow['Solver']} (N={lRow['N']})" if lRow['Solver'] == 'DNS' else f"{lRow['Solver']} (N={lRow['N']}, $C_K$={lRow['Ck']})"
            lRow['label'] = sLab; data.append(lRow)
    
    # Dataframe
    dfData = pd.DataFrame(data); # dfData = dfData.sort_values(by=['Solver', 'N']).reset_index(drop=True)

    # Styles
    plt.figure(figsize=(9, 6))
    stylePlt = {
        'DNS': {'marker': 'x', 'linestyle': '-'},
        'LES': {'marker': 'o', 'linestyle': '--'},
        'DEFAULT': {'marker': 'x', 'linestyle': '-'},
    }

    # Plotting
    for i, row in dfData.iterrows():
        # Data
        vVal = int(row['N']); xVec = np.arange(1, vVal+1); 
        yCol = [str(k) for k in range(1, vVal+1)]; yVec = row[yCol].values.astype(float)

        # Plot
        style = stylePlt.get(row['Solver'], stylePlt['DEFAULT']) 
        plt.plot(xVec, yVec, label=row['label'], marker=style['marker'], linestyle=style['linestyle'], ms=3, linewidth=1.5)

    # Slope
    xVec = np.arange(1, 100+1); yVec = [1/int(k)/int(k) for k in xVec]
    plt.plot(xVec, yVec, label='Slope m=2', linestyle='-.', linewidth=1.5)

    # Configuration
    plt.xscale('log'); plt.yscale('log')
    plt.xticks(fontsize=PG.iFontSizeTicks); plt.yticks(fontsize=PG.iFontSizeTicks)
    plt.xlabel('$k$', fontsize=14); plt.ylabel('$E_k$', fontsize=PG.iFontSizeLabel)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(); plt.tight_layout()

    # Save Plot
    writePath = pathPlt / "Burgers"; writePath.mkdir(parents=True, exist_ok=True)
    tempName = "Plot_Reference.png"; 
    if not os.path.exists(writePath / tempName):
        print("Exporting image ..."); plt.savefig(writePath / tempName); print(f"File saved to: {writePath / tempName}")

    return


def plotNumericalBurgers(fName:str, varName:str): # 1D Burgers Numerical Study
    # Control
    data = []

    # Esto tiene que detectar el archivo y generar los Plots
    # Plot_0 = Energy cascade -- Plot_1 = Velocity in Fourier Domain

    # Import Data
    pDir = [d for d in pathRes.glob(f"*{fName}*") if d.is_dir()]
    for vDir in pDir:
        # Find file
        pFiles = list(vDir.glob("*.csv"))
        for pFile in pFiles:
            # Import File
            df = pd.read_csv(pFile)
            if df.empty: continue
            lRow = df.iloc[-1].to_dict()

            # Data
            argVar = re.search(rf'{varName}(\d+(?:\.\d+)?)', vDir.name); lRow[varName] = float(argVar.group(1))
            argSolver = re.search(r'(DNS|LES)', vDir.name); lRow['Solver'] = argSolver.group(1)
            lRow['label'] = f"{lRow['Solver']} ({varName}={lRow[varName]})"
            data.append(lRow)

    # Dataframe
    dfData = pd.DataFrame(data); dfData = dfData.sort_values(by=['Solver', varName]).reset_index(drop=True)

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
            val['twin'].set_ylabel(f"$E_k$ (Diverged)", fontsize=10)
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
        # val['plot'].set_xticks(fontsize=10); val['plot'].set_yticks(fontsize=10)
        val['plot'].tick_params(axis='both', labelsize=10)
        val['plot'].set_xlabel('$k$', fontsize=14); val['plot'].set_ylabel(f"$E_k$", fontsize=14)
        val['plot'].grid(True, which="both", ls="--", alpha=0.5)
        plt.tight_layout()

        # Legend
        handPlot, labPlot = val['plot'].get_legend_handles_labels()
        val['plot'].legend(handPlot + handTwin, labPlot + labTwin, loc='lower left', fontsize='small');

    # Save Plot
    writePath = pathPlt / "Burgers"; writePath.mkdir(parents=True, exist_ok=True)
    tempName = "Plot_" + fName + ".png"; 
    if not os.path.exists(writePath / tempName):
        print("Exporting image ..."); plt.savefig(writePath / tempName); print(f"File saved to: {writePath / tempName}")

    return
