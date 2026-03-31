import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animate


# Data
filePath = r"C:\Users\gonce\Documents\Master - UPC\0. TFM\HMT\TestData\20260330122723_data_implicit\Probe_1_Map.csv"
data = pd.read_csv(filePath)
vTitle = data.columns.values

aTemp, bTemp = [], []
k, tempVal, fps = 0, 0, 40

# Read Line
for i in range(1, len(vTitle)-1):

    # Coordinates
    tVec = vTitle[i].split(" "); tVec = [float(x) for x in tVec]

    # Check
    if tVec[0] == tempVal:
        bTemp.append(data.iloc[0][i])
    else:
        bTemp = [float(x) for x in bTemp]
        aTemp.append(bTemp)
        bTemp = [data.iloc[0][i]]
    
    # Control
    tempVal = tVec[0]

# Single Plot
# plt.figure(1); plt.imshow(aTemp, cmap='jet')

# Start Plot
fig, ax = plt.subplots()
im = ax.imshow(aTemp, cmap='jet', interpolation='bilinear')
fig.colorbar(im, label="Temperature (°C)")
ax.set_title("Temperature Evolution: Frame 0")

def update(frame):

    # nonlocal k

    # k += 1

    aTemp, bTemp = [], []
    tempVal = 0

    # Read Line
    for i in range(1, len(vTitle)-1):

        # Coordinates
        tVec = vTitle[i].split(" "); tVec = [float(x) for x in tVec]

        # Check
        if tVec[0] == tempVal:
            bTemp.append(data.iloc[int(frame)][i])
        else:
            bTemp = [float(x) for x in bTemp]
            aTemp.append(bTemp)
            bTemp = [data.iloc[int(frame)][i]]
        
        # Control
        tempVal = tVec[0]
    # aTemp.append(bTemp)


    # Update
    aTemp = np.array(aTemp)
    im.set_array(aTemp)
    ax.set_title(f"Temperature Evolution: Frame {frame}")

    return [im]


# Animation
ani = animate.FuncAnimation(fig, update, frames=len(data.index), interval=1, blit=True)

# End
plt.show()
