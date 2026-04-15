import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animate

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
        vTime.append(row[0])

        # Temperature
        for i in range(N[0]):
            for j in range(N[1]):
                
                # Save Node
                k = i * N[1] + j; bTemp.append(row[k+1])
            
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

    # Animation
    ani = animate.FuncAnimation(fig, update, frames=len(frames), interval=0.5, blit=True, repeat=False)

    # Save Video
    if not os.path.exists(filePath + "\\Animation_" + str(idAnimation+1) + ".mp4"):
        print("Exporting video ...")
        idAnimation += 1
        ani.save(filePath + "\\Animation_" + str(idAnimation) + ".mp4", writer='ffmpeg', fps=30)
        print("File saved to: " + filePath + "\\Animation_" + str(idAnimation) + ".mp4")

def createSnapshot(filePath:str, frames:list, vTime:list, tStep:float = -1):

    global idPlot

    # Index
    if tStep == -1: tPos = int(tStep)
    else:
        # Find index and plot that step
        tPos = min(range(len(vTime)), key=lambda i: abs(vTime[i] - tStep))

    # Plot Map
    plt.figure(); im = plt.imshow(frames[tPos], origin='lower')
    plt.colorbar(im, label="Temperature (°C)")
    plt.xlabel('y coordinate'); plt.ylabel('x coordinate')
    plt.title(f'Temperature map @ t = {tStep}')
    
    # Isotherms
    contours = plt.contour(frames[tPos], colors='black', origin='lower', linewidths=0.8, alpha=0.7)
    plt.clabel(contours, inline=True, fontsize=8, fmt='%.1f °C')

    # Save Plot
    if not os.path.exists(filePath + "\\Plot_" + str(idPlot+1) + ".png"):
        print("Exporting image ...")
        idPlot += 1
        plt.savefig(filePath + "\\Plot_" + str(idPlot) + ".png")
        print("File saved to: " + filePath + "\\Plot_" + str(idPlot) + ".png")

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
    if not os.path.exists(filePath + "\\Plot_" + str(idPlot + 1) + ".png"):
        print("Exporting image ...")
        idPlot += 1
        plt.savefig(filePath + "\\Plot_" + str(idPlot) + ".png")
        print("File saved to: " + filePath + "\\Plot_" + str(idPlot) + ".png")

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
    if not os.path.exists(filePath + "\\Plot_" + str(idPlot + 1) + ".png"):
        print("Exporting image ...")
        idPlot += 1
        plt.savefig(filePath + "\\Plot_" + str(idPlot) + ".png")
        print("File saved to: " + filePath + "\\Plot_" + str(idPlot) + ".png")


########## Types of Plots ##########
### Animation (Animate Map): Reads coordinates from file header
### Snapshot (Plot Map): If timeStep not specified, assumes last frame
### Profile (Plot Temperature Profile): If timeStep not specified assumes last frame, needs separate file with a single row of temperatures
### Point (Plot Point Evolution): Plot data from Plot_0_Point.

##### Directory #####
filePath = r"C:\Users\gonce\Documents\Master - UPC\0. TFM\HMT\TestData\20260415083501_test4_cranknicolson" # Change path here

##### Plotting #####

# # Parse Data
# fileName = "\\Probe_1_Map.csv"
# frames, vTime, xVec, yVec = getFrames(filePath + fileName)

# Plot
# createAnimation(filePath, frames, vTime)
# createSnapshot(filePath, frames, vTime, 2000)
# createSnapshot(filePath, frames, vTime, 3000)
# createSnapshot(filePath, frames, vTime, 4000)
# createSnapshot(filePath, frames, vTime, 5000)


# Parse Data
fileName = "\\Probe_2_Map.csv"
frames, vTime, xVec, yVec = getFrames(filePath + fileName)

# Plot
createProfile(filePath, xVec, frames, vTime, True, 3)
# createProfile2(filePath, yVec, frames, vTime, True, 3)
# createPoint(filePath, fileName, 1000)


# Show Plots (Comment if undesired)
plt.show()
