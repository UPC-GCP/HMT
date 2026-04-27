import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animate
from pathlib import Path

def getFrames(fPath:str):

    tStart = time.time()

    # Read Data
    data = pd.read_csv(fPath)
    vTitle = data.columns.values

    # Control
    vRet, vTime = [], []
    N, aTemp, bTemp = [], [], []

    # Dimensions
    i = 1; xPos = vTitle[i].split(' ')[0]
    while vTitle[i+1].split(' ')[0] == xPos: i += 1
    N = [int((len(vTitle)-1) / i), i]
    
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

    return vRet, vTime


########## Plot Map ##########
# Parse Data
basePath = Path.cwd()
dirPath = basePath / "TestData" / "20260422024433_4Materials_implicit"
fileName = dirPath / "Probe_1_Map.csv"
frames, vTime = getFrames(fileName)

# ##### Single Plot #####
# plt.figure(1); plt.imshow(frames[0], cmap='jet')
# plt.figure(2); plt.imshow(frames[100], cmap='jet')
# plt.figure(3); plt.imshow(frames[200], cmap='jet')
# plt.figure(4); plt.imshow(frames[300], cmap='jet')
# plt.figure(5); plt.imshow(frames[400], cmap='jet')
# plt.figure(6); plt.imshow(frames[500], cmap='jet')

fig, ax = plt.subplots()
im = ax.imshow(frames[0], cmap='jet', interpolation='bilinear')
cb = fig.colorbar(im, label="Temperature (°C)")
ax.set_title(f"Temperature Evolution: Time {vTime[0]:.2f} s")

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
fileName = dirPath / "Animation_1.mp4"
if not os.path.exists(fileName):
    print("Exporting video ...")
    ani.save(fileName, writer='ffmpeg', fps=30)
    print(f"File saved to: {fileName}")

