### Imports
import os
import re
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

### Self-Imports
import plotImport as PI

##### Directory #####
pathRes = Path.cwd() / 'ioRes'
pathPlt = Path.cwd() / 'ioPlots'

##### Plotting #####
def plotNumericalBurgers(fName:str, varName:str): # 1D Burgers Numerical Study
    
    # Control
    data = []

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
    print(dfData)

    # Styles
    fig, (axDNS, axLES) = plt.subplots(1, 2, figsize=(18, 6), sharey=True)
    axDNS_t = axDNS.twinx(); axLES_t = axLES.twinx();
    # axDNS_t.set_xscale('log'); axDNS_t.set_yscale('log')
    # axLES_t.set_xscale('log'); axLES_t.set_yscale('log')
    stylePlt = {
        'DNS': {'plot': axDNS, 'twin': axDNS_t, 'bTwin': False, 'marker': 'x', 'linestyle': '--'},
        'LES': {'plot': axLES, 'twin': axLES_t, 'bTwin': False, 'marker': 'x', 'linestyle': '--'},
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
        style = stylePlt.get(row['Solver'], stylePlt['DEFAULT'])
        if np.max(yVec) > 1:
            tempax = style['twin']; templinestyle = ':'
            templabel = f"{row['label']} (Diverged)"
            stylePlt[row['Solver']]['bTwin'] = True
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
            val['twin'].set_ylabel(f"$E_k$ (Diverged)")
            handTwin, labTwin = val['twin'].get_legend_handles_labels()
        else:
            val['twin'].yaxis.set_visible(False)

        # Configuration (Plot)
        val['plot'].set_xscale('log'); val['plot'].set_yscale('log');
        val['plot'].set_xlabel(varName); val['plot'].set_ylabel(f"$E_k$")
        val['plot'].grid(True, which="both", ls="--", alpha=0.5)

        # Legend
        handPlot, labPlot = val['plot'].get_legend_handles_labels()
        val['plot'].legend(handPlot + handTwin, labPlot + labTwin, loc='lower left', fontsize='small');


        # # Slope
        # if True:
        #     varVal = (100 if key == 'DNS' else 30)
        #     xVec = np.arange(1, varVal+1); yVec = [1/int(k) for k in xVec]
        #     val['plot'].plot(xVec, yVec, label='Slope m=2', linestyle='-', linewidth=1.5)

    # Save Plot
    writePath = pathPlt / "Burgers"; writePath.mkdir(parents=True, exist_ok=True)
    tempName = "Plot_" + fName + ".png"; 
    if not os.path.exists(writePath / tempName):
        print("Exporting image ..."); plt.savefig(writePath / tempName); print(f"File saved to: {writePath / tempName}")

    return

##### Run #####
plotNumericalBurgers(sys.argv[1], sys.argv[2])
