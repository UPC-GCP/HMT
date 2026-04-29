import sys
import numpy as np
import PlotLib as PL
from pathlib import Path

########## Main Plotter File ##########
### Use this file to generate the desired plots.
### Plot functions should be added in PlotLib.py

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
if len(sys.argv) != 5: print("Arguments passed incorrectly. (Expected values: Path, bPoint, bAnimate, vSnapshot)"); quit()


### Plot Points
if bool(int(sys.argv[2])): PL.createPoint(dirPath)

### Plot Maps
if bool(int(sys.argv[3])) or len(sys.argv[4]) > 2:
    # Parse Data
    fileName = dirPath / "Probe_1_Map.csv"
    frames, vTime, xVec, yVec = PL.getFrames(fileName)
    cleanFrames = np.array([[x[1:-1] for x in y[1:-1]] for y in frames])

    ### Animate Map
    if bool(int(sys.argv[3])): PL.createAnimation(dirPath, cleanFrames, vTime)
    
    ### Plot Snapshots
    if len(sys.argv[4]) > 2:
        for iSnap in sys.argv[4][1:-1].split(","): PL.createSnapshot(dirPath, cleanFrames, vTime, float(iSnap))


