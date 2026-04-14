import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def getPlots(fPath:str):

    # Get all vectors so you can plot them with coordinates as label

    tStart = time.time()

    # Read Data
    data = pd.read_csv(fPath)
    vTitle = data.columns.values

    # Control
    vRet, vTime = [], data['Time']

    for i in range(len(vTitle)):
        if i == 0: continue
        vRet.append(data[vTitle[i]])

    return vRet, vTime, vTitle[1:]


########## Plot Temperature Evolution ##########
# Parse Data
filePath = r"C:\Users\gonce\Documents\Master - UPC\0. TFM\HMT\TestData\20260408055922_data_crank-nicolson"
fileName = "\\Probe_3_Bug.csv"
xPlot, xTime, xLabel = getPlots(filePath + fileName)

# Plot
fig = plt.figure(1); fig.set_figheight(4); fig.set_figwidth(8)
for i in range(len(xPlot)):
    plt.plot(xTime, xPlot[i], label=xLabel[i])
plt.grid(which="both", alpha=0.2); plt.minorticks_on();
plt.legend(); plt.xlabel('Time (s)'); plt.ylabel('Temperature (°C)')

# Save Plot
if not os.path.exists(filePath + "\\Plot_1.png"):
    print("Exporting image ...")
    plt.savefig(filePath + "\\Plot_1.png")
    print("File saved to: " + filePath + "\\Plot_1.png")

# End
plt.show()
