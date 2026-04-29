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

def plotAnalytic(xVec:list, case:int):
    
    if case == 0: # Dirichlet - Dirichlet

        # Test Variables
        Tl, Tr = 20, 100
        qV, V = 300000, 1
        lamb = 400

        # Coefficients
        C2 = Tl; C1 = 0.5 * qV * V / lamb + Tr - Tl

        # Return
        yRet = [float(C2 + C1*x - 0.5 * qV * V * x * x / lamb) for x in xVec]

    elif case == 1: # Neumann - Dirichlet
        
        # Test Variables
        qNeu, TDir = 10000, 20
        qV, V = 300000, 1
        lamb = 400

        # Coefficients
        C1 = - qNeu/lamb
        C2 = TDir + 0.5 * qV / lamb - C1

        # Return
        yRet = [float(C2 + C1*x - 0.5 * qV * x * x / lamb) for x in xVec]

    elif case == 2: # Convection - Dirichlet

        # Test Variables
        Tg, alpha = 80, 2000; TDir = 20
        qV, V = 300000, 1
        lamb = 400

        # Coefficients
        C2 = (TDir + 0.5*qV/lamb + alpha*Tg/lamb) / (1 + alpha/lamb)
        C1 = - alpha * (Tg - C2) / lamb

        # Return
        yRet = [float(C2 + C1*x - 0.5 * qV * V * x * x / lamb) for x in xVec]
    
    elif case == 3: # Dirichlet - Neumann
        
        # Test Variables
        qNeu, TDir = 15000, 100
        qV, V = 300000, 1
        lamb = 400

        # Coefficients
        C2 = TDir
        C1 = qV / lamb - qNeu / lamb

        # Return
        yRet = [float(C2 + C1*x - 0.5 * qV * V * x * x / lamb) for x in xVec]
    
    elif case == 4: # Dirichlet - Convection
        
        # Test Variables
        Tg, alpha = 10, 4000; TDir = 100
        qV, V = 300000, 1
        lamb = 400

        # Coefficients
        C2 = TDir
        C1 = (alpha*Tg/lamb + (1 + 0.5*alpha/lamb)*(qV/lamb) - alpha*C2/lamb) / (1 + alpha/lamb)

        # Return
        yRet = [float(C2 + C1*x - 0.5 * qV * V * x * x / lamb) for x in xVec]
    
    return yRet

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

def createAnimation(filePath:str, frames:list, vTime:list):

    global idAnimation

    # Plot Map
    fig, ax = plt.subplots()
    im = ax.imshow(frames[0], cmap='jet', interpolation='bilinear')
    cb = fig.colorbar(im, label="Temperature (°C)")
    ax.set_title(f"Temperature Evolution: Time {vTime[0]:.2f} s")

    # Update Function
    def update(frame):
    
        # Control
        arrData = frames[frame]

        # Update Plot
        im.set_array(arrData); im.set_clim(np.min(arrData), np.max(arrData))
        cb.update_normal(im)
        ax.set_title(f"Temperature Evolution: Time {vTime[frame]:.2f} s");

        fig.canvas.draw()

        return [im]
    
    # Control
    tStart = time.time()
    print(f"Parsing file: {str(filePath).split("/")[-1]}")


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

def createSnapshot(filePath:str, frames:list, vTime:list, tStep:float = -1):

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

    # Plot Map
    plt.figure(); im = plt.imshow(frames[tPos], origin='lower')
    plt.colorbar(im, label="Temperature (°C)")
    plt.xlabel('y coordinate'); plt.ylabel('x coordinate')
    plt.title(f'Temperature map @ t = {tStep}')
    
    # Isotherms
    contours = plt.contour(frames[tPos], levels=levels, colors='black', origin='lower', linewidths=0.8, alpha=0.7)
    plt.clabel(contours, inline=True, fontsize=8, fmt='%.1f °C')

    # Save Plot
    tempName = "Plot_" + str(idPlot + 1) + ".png"
    if not os.path.exists(filePath / tempName):
        print("Exporting image ...")
        idPlot += 1
        plt.savefig(filePath / tempName)
        print(f"File saved to: {filePath / tempName}")

def createProfile(filePath:str, xVec:list, frames:list, vTime:list, bAnal=False, iAnal=0, tStep:float = -1):

    global idPlot

    # Index
    if tStep == -1: tPos = -1
    else:
        # Find index and plot that step
        tPos = min(range(len(vTime)), key=lambda i: abs(vTime[i] - tStep))

    # Plot Map
    plt.figure()
    plt.plot(xVec, frames[tPos], 'r',label='Simulation Results')
    if bAnal:
        aSol = plotAnalytic(xVec, iAnal)
        plt.plot(xVec, aSol, 'b', label='Analytical Solution')
    plt.xlabel('Position (m)'); plt.ylabel('Temperature (°C)')
    plt.legend(); plt.grid(which="both", alpha=0.2); plt.minorticks_on()

    # Save Plot
    tempName = "Plot_" + str(idPlot + 1) + ".png"
    if not os.path.exists(filePath / tempName):
        print("Exporting image ...")
        idPlot += 1
        plt.savefig(filePath / tempName)
        print(f"File saved to:  {filePath / tempName}")

def createProfile2(filePath:str, xVec:list, frames:list, vTime:list, bAnal=False, iAnal=0, tStep:float = -1):

    global idPlot

    # Index
    if tStep == -1: tPos = -1
    else:
        # Find index and plot that step
        tPos = min(range(len(vTime)), key=lambda i: abs(vTime[i] - tStep))

    # Clean Vector
    for nFrame in frames[tPos]:
        newFrames = [float(x) for x in nFrame]

    # Plot Map
    plt.figure()
    plt.plot(xVec, newFrames, 'r',label='Simulation Results')
    if bAnal:
        aSol = plotAnalytic(xVec, iAnal)
        plt.plot(xVec, aSol, 'b', label='Analytical Solution')
    plt.xlabel('Position (m)'); plt.ylabel('Temperature (°C)')
    plt.legend(); plt.grid(which="both", alpha=0.2); plt.minorticks_on()

    # Save Plot
    tempName = "Plot_" + str(idPlot + 1) + ".png"
    if not os.path.exists(filePath / tempName):
        print("Exporting image ...")
        idPlot += 1
        plt.savefig(filePath / tempName)
        print(f"File saved to: {filePath / tempName}")

def createPoint(filePath:str, xVec:list, vTime:list):

    pass


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
### iAnalytical: Create comparison with analytical Solution


##### Directory #####
dirPath = Path.cwd() / "TestData" / sys.argv[1]
if len(sys.argv) != 7: print("Arguments passed incorrectly. (Expected values: Path, bPoint, bAnimate, vSnapshot, iAnalytical, iAxis)"); quit()


### Plot Points
if bool(int(sys.argv[2])):
    print("Plot Probe_Point data")

    # Parse Data
    # fileName = dirPath / "Probe_0_Point.csv"
    # frames, vTime, xVec, yVec = getFrames(fileName)

    # Plot Data
    # createPoint()


if bool(int(sys.argv[3])) or len(sys.argv[4]) > 2:
    
    # Parse Data
    fileName = dirPath / "Probe_1_Map.csv"
    frames, vTime, xVec, yVec = getFrames(fileName)

    ### Animate Map   
    if bool(int(sys.argv[3])): 

        # Plot
        createAnimation(dirPath, frames, vTime)
    
    ### Plot Snapshots
    for iSnap in sys.argv[4][1:-1].split(","):
        
        # Filter
        if len(iSnap) == 0: break

        # Plot
        createSnapshot(dirPath, frames, vTime, float(iSnap))


### Plot Analytical
if int(sys.argv[5]) != -1:

    print("Plot analytical solution comparison")

    # Parse Data
    fileName = dirPath / "Probe_2_Map.csv"
    frames, vTime, xVec, yVec = getFrames(fileName)

    # Plot
    if int(sys.argv[6]) == 0: createProfile(dirPath, xVec, frames, vTime, True, int(sys.argv[5]))
    elif int(sys.argv[6]) == 1: createProfile2(dirPath, yVec, frames, vTime, True, int(sys.argv[5]))

quit()

##### Plotting #####

# Parse Data
fileName = dirPath / "Probe_1_Map.csv"
frames, vTime, xVec, yVec = getFrames(fileName)

# Plot Maps
if int(sys.argv[3]) == 1: createAnimation(dirPath, frames, vTime)
createSnapshot(dirPath, frames, vTime, 2000)
createSnapshot(dirPath, frames, vTime, 3000)
createSnapshot(dirPath, frames, vTime, 4000)
createSnapshot(dirPath, frames, vTime, 5000)


# Parse Data
fileName = dirPath / "Probe_2_Map.csv"
frames, vTime, xVec, yVec = getFrames(fileName)

# Plot Lines
# if int(sys.argv[2]) != -1: createProfile(dirPath, xVec, frames, vTime, True, int(sys.argv[2])) # Only fixed one
if int(sys.argv[2]) != -1: createProfile2(dirPath, yVec, frames, vTime, True, int(sys.argv[2]))

# createPoint(filePath, fileName, 1000)

