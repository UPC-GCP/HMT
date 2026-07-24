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
### bVelocity: Create velocity field map and snapshot
### bNumerical: Plot numerical study
### sVar: Variable for numerical study


# Safety Check
# if len(sys.argv) != 5 and len(sys.argv) != 6 and len(sys.argv) != 8: print("Arguments passed incorrectly. (Expected values: Path, bPoint, bAnimate, vSnapshot)"); print("Additional arguments available. (Optional values: bVelocity, iNumerical, sVar)"); quit()
# if len(sys.argv) != 0: print("Arguments passed incorrectly.")


##### Directory #####
# dirPath = Path.cwd() / "ioRes" / sys.argv[1]
dirPath  = Path.cwd() / 'ioRes'
plotPath = Path.cwd() / 'ioPlots'


# ### Plot Points
# if bool(int(sys.argv[2])): PL.createPoint(dirPath)

# ### Plot Maps
# if bool(int(sys.argv[3])) or len(sys.argv[4]) > 2:
#     # Parse Data
#     fileName = dirPath / "Probe_1_Map.csv"
#     frames, vTime, xVec, yVec = PL.getFrames(fileName)
#     cleanFrames = np.array([[x[1:-1] for x in y[1:-1]] for y in frames])
#     # cleanFrames = frames

#     ### Animate Map
#     if bool(int(sys.argv[3])): PL.createAnimation(dirPath, cleanFrames, vTime, outPath=plotPath / sys.argv[1])

#     ### Plot Snapshots
#     if len(sys.argv[4]) > 2:
#         for iSnap in sys.argv[4][1:-1].split(","): PL.createSnapshot(dirPath, cleanFrames, vTime, float(iSnap), outPath=plotPath / sys.argv[1])

# if len(sys.argv) == 5: quit()

### Convection Study
# PL.numStudyConv(dirPath)

### Navier-Stokes Study
PL.numStudyNS(dirPath)

# ### Plot NS Map
# if bool(int(sys.argv[5])):
#     # Read velocity data and create new vector for total magnitude
#     # Use createSnapshot() and pass the new velocity map
#     pass

# if len(sys.argv) == 6: quit()

# ### Numerical Studies
# if int(sys.argv[6]) == 0: PL.numStudyDiff(fileName, sys.argv[7])
# elif int(sys.argv[6]) == 1: PL.numStudyConv(fileName, sys.argv[7])
# elif int (sys.argv[6]) == 2: PL.numStudyNS(fileName, sys.argv[7])
# else: print("Study selection not recognized...")
    

