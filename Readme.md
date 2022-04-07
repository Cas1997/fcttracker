# ALICE 3 Forward Tracker

### 1. Simulation

`$ o2-sim -m FCT -e TGeant3 -g pythia8 -n 1000`

Output: `o2sim_HitsFCT.root`

### 2. Tracking
`$ root.exe -q fctTracker.C+`

Output: `fcttracks.root`

### 3. Assessment histograms
#### Layer hit densities

`$ root.exe -b -q fctOccupancy.C+`

Output: `fctOccupancy.root`

#### Tracking evalutaion
`$ root.exe -q -b FCTTrackerChecker.C+`

Output: `Fittercheck_fcttracks.root`
