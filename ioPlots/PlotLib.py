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

def getLastFrame(fPath):
    """Read header and last data row from a probe CSV line-by-line (O(1) memory)."""
    with open(fPath, 'rb') as f:
        header = f.readline().decode().rstrip()
        lastLine = b''
        for line in f:
            if line.strip(): lastLine = line
    vTitle = header.split(',')

    # If the Field probe header is missing its trailing newline, the first data row
    # is merged onto the header line. Trim vTitle to only coordinate columns
    # (format "X Y" — contain a space); plain data values never do.
    lastCoord = max((idx for idx in range(1, len(vTitle)) if ' ' in vTitle[idx]), default=0)
    vTitle = vTitle[:lastCoord + 1]

    vals   = lastLine.decode().rstrip().split(',')
    t      = float(vals[0])
    data   = np.array([float(v) for v in vals[1:] if v])

    i = 1; xFirst = vTitle[i].split(' ')[0]
    xRet = [float(vTitle[i].split(' ')[0])]; yRet = [float(vTitle[i].split(' ')[1])]
    while i + 1 < len(vTitle) and vTitle[i+1].split(' ')[0] == xFirst:
        i += 1; yRet.append(float(vTitle[i].split(' ')[1]))
    nY = i; nX = (len(vTitle) - 1) // nY
    for j in range(1, nX): xRet.append(float(vTitle[1 + j*nY].split(' ')[0]))

    return data.reshape(nX, nY), t, xRet, yRet

def createAnimation(filePath:str, frames:list, vTime:list, outPath:str=None): # TRANSPOSE PENDING

    global idAnimation

    # Plot Map
    fig, ax = plt.subplots()
    im = ax.imshow(np.transpose(frames[0]), cmap='jet', interpolation='bilinear', origin='lower', extent=[-1, 1, 0, 1])
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
    writePath = Path(outPath) if outPath else Path(filePath)
    writePath.mkdir(parents=True, exist_ok=True)
    tempName = "Animation_" + str(idAnimation+1) + ".mp4"
    if not os.path.exists(writePath / tempName):
        print("Exporting video ...")
        idAnimation += 1
        ani.save(writePath / tempName, writer='ffmpeg', fps=30)
        print(f"File saved to: {writePath / tempName}")

    print(f"Elapsed Time: {time.time() - tStart:.3f}")

def createSnapshot(filePath:str, frames:list, vTime:list, tStep:float = -1, extent=None, label="Value", outPath:str=None): # FUNCTIONAL

    global idPlot

    # Index
    if tStep == -1: tPos = -1
    else:
        tPos = min(range(len(vTime)), key=lambda i: abs(vTime[i] - tStep))

    # Levels — predefined for known diffusion timesteps, auto otherwise
    levelsMap = {2000: [15,16,17,18,19,20,21,22], 3000: [18.5,19,19.5,20,20.5,21,21.5,22,22.5], 4000: [21,22,23,24,25,26,27], 5000: [23,24,25,26,27,28,29,30,31,32]}
    levels = levelsMap.get(tStep, None)

    # Frame
    framePlot = np.transpose(frames[tPos])

    # Extent
    if extent is None: extent = [0, 1.1, 0, 0.8]

    # Plot Map
    plt.figure(); im = plt.imshow(framePlot, cmap='jet', extent=extent, origin='lower', aspect='auto')
    plt.colorbar(im, label=label)
    plt.xlabel('x (m)'); plt.ylabel('y (m)')
    plt.title(f'{label} map @ t = {vTime[tPos]:.4f} s')

    # Isotherms — only when levels are defined
    if levels is not None:
        contours = plt.contour(framePlot, levels=levels, colors='black', origin='lower', linewidths=0.8, alpha=0.7, extent=extent)
        plt.clabel(contours, inline=True, fontsize=8, fmt='%.1f')

    # Save Plot
    writePath = Path(outPath) if outPath else Path(filePath)
    writePath.mkdir(parents=True, exist_ok=True)
    tempName = "Plot_" + str(idPlot + 1) + ".png"
    if not os.path.exists(writePath / tempName):
        print("Exporting image ...")
        idPlot += 1
        plt.savefig(writePath / tempName)
        print(f"File saved to: {writePath / tempName}")


def createPoint(filePath:str, outPath:str=None): # FUNCTIONAL

    global idPlot
    writePath = Path(outPath) if outPath else Path(filePath)
    writePath.mkdir(parents=True, exist_ok=True)

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
    if not os.path.exists(writePath / tempName):
        print("Exporting image ..."); idPlot += 1
        plt.savefig(writePath / tempName)
        print(f"File saved to: {writePath / tempName}")

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
    if not os.path.exists(writePath / tempName):
        print("Exporting image ..."); idPlot += 1
        plt.savefig(writePath / tempName)
        print(f"File saved to: {writePath / tempName}")

    
def numStudyDiff(fileName:str, sVar:str) :

    # Range
    if sVar == 'dt':
        rCase = [1,2,3,4,9,10,11,12,17,18,19,20]
        rVar = [0.5, 1, 2, 5]
    elif sVar == 'N':
        rCase = [5,6,7,8,13,14,15,16,21,22,23,24]
        rVar = [50, 100, 200, 400]
    else: print('Variable not recognized ...'); quit
    
    # Directory
    dataPath = Path.cwd() / "TestData"
    outPath  = Path(__file__).parent / "DiffStudy"
    outPath.mkdir(parents=True, exist_ok=True)
    vRes, vIter = [], []
    fileList = [f for f in sorted(os.listdir(dataPath)) if "Case" in f]

    # File Loop
    for file in fileList:

        # Filter
        if not int(file[20:22]) in rCase: continue

        # Read Data
        data = pd.read_csv(dataPath / file / fileName)

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
        if not os.path.exists(outPath / tempName):
            print("Exporting image ...")
            plt.savefig(outPath / tempName)
            print(f"File saved to: {outPath / tempName}")
            print(f"File saved to: {dirPath / tempName}")


def numStudyConv(dirPath: Path):

    # Reference data: Case 0: rho/Gamma=10, Case 1: rho/Gamma=1000, Case 2: rho/Gamma=1e6
    # xAnal covers the Neumann half of the South wall (x in [0, 1])
    xAnal  = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    yAnal  = [
        [1.989, 1.402, 1.146, 0.946, 0.775, 0.621, 0.480, 0.349, 0.227, 0.111, 0.000],
        [2.0000, 1.9990, 1.9997, 1.9850, 1.8410, 0.9510, 0.1540, 0.0010, 0.0000, 0.0000, 0.0000],
        [2.000, 2.000, 2.000, 1.999, 1.964, 1.000, 0.036, 0.001, 0.000, 0.000, 0.000],
    ]
    rhoGamma = [10, 1000, 1000000]
    nCases   = len(rhoGamma)
    schemes  = ['UDS', 'CDS']
    aColor   = {'UDS': 'b', 'CDS': 'r'}

    # Output folder for comparison plots
    outPath = Path(__file__).parent / "ConvStudy"
    os.makedirs(outPath, exist_ok=True)

    # Build case→directory map: caseDirs[scheme][sCase] → Path
    # UDS: exConv{1,2,3}; CDS: exConv{4,5,6}
    caseDirs = {s: {} for s in schemes}
    for sCase in range(nCases):
        for scheme, num in [('UDS', sCase+1), ('CDS', sCase+4)]:
            matches = sorted([d for d in os.listdir(dirPath) if f"exConv{num}_" in d])
            if matches: caseDirs[scheme][sCase] = dirPath / matches[-1]

    # ---- 1. Colormap snapshot (all 6 cases, saved per case directory) ----
    for scheme in schemes:
        for sCase in range(nCases):
            if sCase not in caseDirs[scheme]: continue
            cDir = caseDirs[scheme][sCase]
            frames, vTime, xVec, yVec = getFrames(cDir / "Probe_1_Map.csv")
            extent = [xVec[0], xVec[-1], yVec[0], yVec[-1]]
            createSnapshot(outPath, frames, vTime, extent=extent, label='Φ')
            plt.close()

    # ---- 2. Mean iterations and residue vs rho/Gamma (one figure each) ----
    vIter = {s: [] for s in schemes}
    vRes  = {s: [] for s in schemes}
    for scheme in schemes:
        for sCase in range(nCases):
            if sCase not in caseDirs[scheme]: continue
            data = pd.read_csv(caseDirs[scheme][sCase] / "Probe_3_Bug.csv")
            vIter[scheme].append(data['lastIter'].mean())
            vRes[scheme].append(data['lastRes'].mean())

    xTicks = [str(r) for r in rhoGamma]
    for metric, vData, yLabel in [('Iterations', vIter, 'Mean Iterations'), ('Residue', vRes, 'Mean Residue')]:
        plt.figure()
        for scheme in schemes:
            if vData[scheme]:
                plt.plot(range(len(vData[scheme])), vData[scheme], color=aColor[scheme], marker='o', label=scheme)
        plt.xticks(range(nCases), xTicks)
        plt.xlabel('ρ/Γ'); plt.ylabel(yLabel)
        plt.legend(); plt.grid(which='both', alpha=0.2); plt.minorticks_on()
        plt.savefig(outPath / f"Conv_{metric}.png"); print(f"Saved: Conv_{metric}.png"); plt.close()

    # ---- 3. South boundary Phi profile vs reference (one figure per material case) ----
    for sCase in range(nCases):
        plt.figure()
        plt.plot(xAnal, yAnal[sCase], 'k--', marker='x', markersize=6, label='Reference')
        for scheme in schemes:
            if sCase not in caseDirs[scheme]: continue
            data = pd.read_csv(caseDirs[scheme][sCase] / "Probe_2_Boundary.csv")
            xCoords = np.array([float(col.split(' ')[0]) for col in data.columns[1:]])
            phiVals = data.iloc[-1, 1:].values.astype(float)
            plt.plot(xCoords, phiVals, color=aColor[scheme], label=scheme)
        plt.xlabel('x (m)'); plt.ylabel('Φ')
        plt.title(f'South Boundary Profile — ρ/Γ = {rhoGamma[sCase]}')
        plt.legend(); plt.grid(which='both', alpha=0.2); plt.minorticks_on()
        plt.savefig(outPath / f"Conv_Boundary_case{sCase}.png"); print(f"Saved: Conv_Boundary_case{sCase}.png"); plt.close()


def numStudyNS(dirPath: Path):

    # Reference data: Ghia et al. (1982), lid-driven cavity
    # u at x=0.5: yU=y-positions (0=bottom, 1=lid), xURef=u-velocity
    yU   = [1, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0]
    xURef = [
        [1, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0],
        [1, 0.75837, 0.68439, 0.61756, 0.55892, 0.29093, 0.16256, 0.02135, -0.11477, -0.17119, -0.32726, -0.24299, -0.14612, -0.10338, -0.09266, -0.08186, 0]
    ]
    # v at y=0.5: xV=x-positions (0=left, 1=right), yVRef=v-velocity
    xV   = [0, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, 0.9688, 1]
    yVRef = [
        [0, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 0.17507, 0.17527, 0.05454, -0.24533, -0.22445, -0.16914, -0.10313, -0.08864, -0.07391, -0.05906, 0],
        [0, 0.18360, 0.19713, 0.20920, 0.22965, 0.28124, 0.30203, 0.30174, 0.05186, -0.38598, -0.44993, -0.23827, -0.22847, -0.19254, -0.15663, -0.12146, 0]
    ]
    Re      = [100, 400]
    schemes = ['UDS', 'CDS']
    aColor  = {'UDS': 'b', 'CDS': 'r'}

    outPath = Path(__file__).parent / "NSStudy"
    os.makedirs(outPath, exist_ok=True)

    # Case dirs: exNS1=UDS/Re100, exNS2=UDS/Re400, exNS3=CDS/Re100, exNS4=CDS/Re400
    caseDirs = {s: {} for s in schemes}
    for sCase in range(2):
        for scheme, num in [('UDS', sCase+1), ('CDS', sCase+3)]:
            matches = sorted([d for d in os.listdir(dirPath) if f"exNS{num}_" in d])
            if matches: caseDirs[scheme][sCase] = dirPath / matches[-1]

    # ---- 1. Velocity magnitude colormaps + vector maps (one per case) ----
    for scheme in schemes:
        for sCase in range(2):
            if sCase not in caseDirs[scheme]: continue
            cDir = caseDirs[scheme][sCase]
            print(f"Reading {scheme} Re={Re[sCase]} ...")
            uFrame, _, uX, uY = getLastFrame(cDir / "Probe_2u_Field.csv")
            vFrame, _, vX, vY = getLastFrame(cDir / "Probe_2v_Field.csv")

            # Interpolate staggered faces to p-grid cell centres:
            # u: pad x=1 wall face (u=0), average adjacent faces in x
            # v: pad y=1 lid face (v=0), average adjacent faces in y
            N     = uFrame.shape[0]
            u_pad = np.vstack([uFrame, np.zeros((1, N))])
            v_pad = np.hstack([vFrame, np.zeros((N, 1))])
            u_p   = 0.5 * (u_pad[:-1] + u_pad[1:])
            v_p   = 0.5 * (v_pad[:, :-1] + v_pad[:, 1:])
            mag   = np.sqrt(u_p**2 + v_p**2)
            pX    = np.array(uX) + 0.5/N   # face positions → cell centres
            pY    = np.array(uY)            # u-mesh y is already at cell centres

            # Velocity magnitude colormap
            plt.figure()
            im = plt.imshow(mag.T, cmap='jet', origin='lower', aspect='equal', extent=[0, 1, 0, 1])
            plt.colorbar(im, label='|V| (m/s)')
            plt.xlabel('x (m)'); plt.ylabel('y (m)')
            plt.title(f'Velocity Magnitude — Re={Re[sCase]}, {scheme}')
            plt.savefig(outPath / f"NS_Magnitude_Re{Re[sCase]}_{scheme}.png")
            plt.close(); print(f"Saved: NS_Magnitude_Re{Re[sCase]}_{scheme}.png")

            # Velocity vector map (subsampled for clarity)
            step = max(1, N // 16)
            XX, YY = np.meshgrid(pX[::step], pY[::step], indexing='ij')
            plt.figure()
            plt.quiver(XX, YY, u_p[::step, ::step], v_p[::step, ::step])
            plt.xlabel('x (m)'); plt.ylabel('y (m)')
            plt.title(f'Velocity Field — Re={Re[sCase]}, {scheme}')
            plt.savefig(outPath / f"NS_Vectors_Re{Re[sCase]}_{scheme}.png")
            plt.close(); print(f"Saved: NS_Vectors_Re{Re[sCase]}_{scheme}.png")

    # ---- 2. Centerline profiles vs Ghia reference (one figure per Re per component) ----
    for sCase in range(2):

        # u along vertical centerline x=0.5
        plt.figure()
        plt.plot(xURef[sCase], yU, 'k--', marker='x', markersize=5, label='Ghia et al.')
        for scheme in schemes:
            if sCase not in caseDirs[scheme]: continue
            uFrame, _, uX, uY = getLastFrame(caseDirs[scheme][sCase] / "Probe_2u_Field.csv")
            iX = np.argmin(np.abs(np.array(uX) - 0.5))
            plt.plot(uFrame[iX, :], uY, color=aColor[scheme], label=scheme)
        plt.xlabel('u (m/s)'); plt.ylabel('y (m)')
        plt.title(f'u at x=0.5 — Re={Re[sCase]}')
        plt.legend(); plt.grid(which='both', alpha=0.2); plt.minorticks_on()
        plt.savefig(outPath / f"NS_u_Re{Re[sCase]}.png"); plt.close()
        print(f"Saved: NS_u_Re{Re[sCase]}.png")

        # v along horizontal centerline y=0.5
        plt.figure()
        plt.plot(xV, yVRef[sCase], 'k--', marker='x', markersize=5, label='Ghia et al.')
        for scheme in schemes:
            if sCase not in caseDirs[scheme]: continue
            vFrame, _, vX, vY = getLastFrame(caseDirs[scheme][sCase] / "Probe_2v_Field.csv")
            jY = np.argmin(np.abs(np.array(vY) - 0.5))
            plt.plot(vX, vFrame[:, jY], color=aColor[scheme], label=scheme)
        plt.xlabel('x (m)'); plt.ylabel('v (m/s)')
        plt.title(f'v at y=0.5 — Re={Re[sCase]}')
        plt.legend(); plt.grid(which='both', alpha=0.2); plt.minorticks_on()
        plt.savefig(outPath / f"NS_v_Re{Re[sCase]}.png"); plt.close()
        print(f"Saved: NS_v_Re{Re[sCase]}.png")

