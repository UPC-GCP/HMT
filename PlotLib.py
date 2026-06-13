import os
import sys
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animate
from pathlib import Path

# Indexing
idAnimation = 0
idPlot = 0

def getFrames(fPath:str):

    tStart = time.time()
    print(f"Parsing file: {str(fPath).split("/")[-1]}")

    # Read Data
    data = pd.read_csv(fPath)
    vTitle = data.columns.values

    # Control
    vRet, vTime, xRet, yRet = [], [], [], []
    N, aTemp, bTemp = [], [], []

    # Dimensions / Coordinates
    i = 1; cordTemp = vTitle[i].split(' ')
    xPos = cordTemp[0]
    xRet.append(cordTemp[0]); yRet.append(cordTemp[1])

    while vTitle[i+1].split(' ')[0] == xPos and i < len(vTitle)-2:
        i += 1
        yRet.append(vTitle[i].split(' ')[1])
    N = [int((len(vTitle)-1) / i), i]
    
    for j in range(1, N[0]):
        xRet.append(vTitle[1 + j*N[1]].split(' ')[0])
    
    xRet = [float(x) for x in xRet]
    yRet = [float(y) for y in yRet]
    
    # Rows Loop
    for index, row in data.iterrows():

        # Time
        vTime.append(row.iloc[0])

        # Temperature
        for i in range(N[0]):
            for j in range(N[1]):
                
                # Save Node
                k = i * N[1] + j; bTemp.append(row.iloc[k+1])
            
            # Save Map
            bTemp = [float(x) for x in bTemp]; aTemp.append(bTemp)
            
            # Control
            bTemp = []

        # Save Frame
        vRet.append(np.array(aTemp))

        # Control
        aTemp = []
    
    print(f"Elapsed Time: {time.time() - tStart:.3f}")

    return vRet, vTime, xRet, yRet

def createAnimation(filePath:str, frames:list, vTime:list): # TRANSPOSE PENDING

    global idAnimation

    # Plot Map
    fig, ax = plt.subplots()
    im = ax.imshow(np.transpose(frames[0]), cmap='jet', interpolation='bilinear', origin='lower', extent=[0, 1.1, 0, 0.8])
    plt.xlabel('Length (m)'); plt.ylabel('Height (m)')
    cb = fig.colorbar(im, label="Temperature (°C)")
    ax.set_title(f"Temperature Evolution: Time {vTime[0]:.2f} s")

    # Update Function
    def update(frame):
    
        # Control
        arrData = np.transpose(frames[frame])

        # Update Plot
        im.set_array(arrData); im.set_clim(np.min(arrData), np.max(arrData))
        cb.update_normal(im)
        ax.set_title(f"Temperature Evolution: Time {vTime[frame]:.2f} s");

        fig.canvas.draw()

        return [im]
    
    # Control
    tStart = time.time()

    # Animation
    ani = animate.FuncAnimation(fig, update, frames=len(frames), interval=0.5, blit=True, repeat=False)

    # Save Video
    tempName = "Animation_" + str(idAnimation+1) + ".mp4"
    if not os.path.exists(filePath / tempName):
        print("Exporting video ...")
        idAnimation += 1
        ani.save(filePath / tempName, writer='ffmpeg', fps=30)
        print(f"File saved to: {filePath / tempName}")

    print(f"Elapsed Time: {time.time() - tStart:.3f}")

def createSnapshot(filePath:str, frames:list, vTime:list, tStep:float = -1): # FUNCTIONAL

    global idPlot

    # Index
    if tStep == -1: tPos = int(tStep)
    else:
        # Find index and plot that step
        tPos = min(range(len(vTime)), key=lambda i: abs(vTime[i] - tStep))

    # Levels
    if tStep == 2000: levels = [15, 16, 17, 18, 19, 20, 21, 22]
    elif tStep == 3000: levels = [18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5]
    elif tStep == 4000: levels = [21, 22, 23, 24, 25, 26, 27]
    elif tStep == 5000: levels = [23, 24, 25, 26, 27, 28, 29, 30, 31, 32]

    # Axis
    framePlot = np.transpose(frames[tPos])

    # Plot Map
    plt.figure(); im = plt.imshow(framePlot, cmap='jet', extent=[0, 1.1, 0, 0.8], origin='lower', aspect='auto')
    plt.colorbar(im, label="Temperature (°C)")
    plt.xlabel('Length (m)'); plt.ylabel('Height (m)')
    plt.title(f'Temperature map @ t = {tStep}')
    
    # Isotherms
    contours = plt.contour(framePlot, levels=levels, colors='black', origin='lower', linewidths=0.8, alpha=0.7, extent=[0, 1.1, 0, 0.8])
    plt.clabel(contours, inline=True, fontsize=8, fmt='%.1f °C')

    # Save Plot
    tempName = "Plot_" + str(idPlot + 1) + ".png"
    if not os.path.exists(filePath / tempName):
        print("Exporting image ...")
        idPlot += 1
        plt.savefig(filePath / tempName)
        print(f"File saved to: {filePath / tempName}")


def createPoint(filePath:str): # FUNCTIONAL
    
    global idPlot

    # Read Data
    print("Parsing file: Probe_0_Point.csv")
    data = pd.read_csv(filePath / "Probe_0_Point.csv")
    vTitle = data.columns.values

    # Benchmark
    TPoint = [0, 1000, 2000, 3000, 4000, 5000]
    DPoint = [[8, 12.03, 16.01, 19.22, 22.02, 24.59], [8, 10.86, 15.26, 19.20, 22.37, 25.52]]
    aColor = ['r', 'b']

    # Plot Figure
    plt.figure()
    for i, sHead in enumerate(vTitle[1:]):
        plt.plot(data['Time'][1:], data[sHead][1:], aColor[i], label=sHead + ' (Model)')
        plt.plot(TPoint, DPoint[i], aColor[i] + '--', label=sHead + ' (Benchmark)')
    plt.xlabel('Time (s)'); plt.ylabel('Temperature (°C)')
    plt.legend(); plt.grid(which='both', alpha=0.2); plt.minorticks_on()

    # Save Plot
    tempName = "Plot_" + str(idPlot+1) + ".png"
    if not os.path.exists(filePath / tempName):
        print("Exporting image ..."); idPlot += 1
        plt.savefig(filePath / tempName)
        print(f"File saved to: {filePath / tempName}")

    # Plot Secondary
    plt.figure()
    for i, sHead in enumerate(vTitle[1:]):
        aMid = len(data[sHead][1:]) // 2
        plt.plot(data['Time'][1:aMid], data[sHead][1:aMid], aColor[i], label=sHead + ' (Model)')
        plt.plot(TPoint, DPoint[i], aColor[i] + '--', label=sHead + ' (Benchmark)')
    plt.xlabel('Time (s)'); plt.ylabel('Temperature (°C)')
    plt.legend(); plt.grid(which='both', alpha=0.2); plt.minorticks_on()

    # Save Secondary
    tempName = "Plot_" + str(idPlot+1) + ".png"
    if not os.path.exists(filePath / tempName):
        print("Exporting image ..."); idPlot += 1
        plt.savefig(filePath / tempName)
        print(f"File saved to: {filePath / tempName}")

    
def createNumericalStudy(fileName:str, sVar:str) :

    # Range
    if sVar == 'dt':
        rCase = [1,2,3,4,9,10,11,12,17,18,19,20]
        rVar = [0.5, 1, 2, 5]
    elif sVar == 'N':
        rCase = [5,6,7,8,13,14,15,16,21,22,23,24]
        rVar = [50, 100, 200, 400]
    else: print('Variable not recognized ...'); quit
    
    # Directory 
    dirPath = Path.cwd() / "TestData"; vRes, vIter = [], []; 
    fileList = [f for f in sorted(os.listdir(dirPath)) if "Case" in f]
    
    # File Loop
    for file in fileList:
        
        # Filter
        if not int(file[20:22]) in rCase: continue

        # Read Data
        data = pd.read_csv(dirPath / file / fileName)
        
        # Log
        vRes.append(data['lastRes'].mean())
        vIter.append(data['lastIter'].mean())

    # Loop
    aColor = ['r', 'b', 'g']; aVal = ['Residue', 'Iterations']
    for i, vPlot in enumerate([vRes, vIter]):
        # Plot
        plt.figure()
        plt.plot(rVar, vPlot[0:4], aColor[0], label='implicit')
        plt.plot(rVar, vPlot[4:8], aColor[1], label='crank-nicolson')
        plt.plot(rVar, vPlot[8:12], aColor[2], label='explicit')
        plt.xlabel(sVar); plt.ylabel(aVal[i]); plt.legend()
        plt.grid(which='both', alpha=0.2); plt.minorticks_on()

        # Save
        tempName = "Plot_" + sVar + "_" + aVal[i] + ".png"
        if not os.path.exists(dirPath / tempName):
            print("Exporting image ...")
            plt.savefig(dirPath / tempName)
            print(f"File saved to: {dirPath / tempName}")


def compConvDiff(sCase:int):
    
    xAnal = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 10]
    if sCase == 0: yAnal = [1.989, 1.402, 1.146, 0.946, 0.775, 0.621, 0.480, 0.349, 0.227, 0.111, 0.000]
    elif sCase == 1: yAnal = [2.0000, 1.9990, 1.9997, 1.9850, 1.8410, 0.9510, 0.1540, 0.0010, 0.0000, 0.0000, 0.0000]
    elif sCase == 2: yAnal = [2.000, 2.000, 2.000, 1.999, 1.964, 1.000, 0.036, 0.001, 0.000, 0.000, 0.000]
