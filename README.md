# Brain Source Localization Project

## Overview

This project involves solving the brain source localization problem using various methods and models. The tasks include generating and visualizing dipole locations, calculating potentials, and applying source localization algorithms.

## 1. Model and Electrode Positioning

### Data
- **File**: `mat.ElecPosXYZ`
- **Description**: Contains normalized electrode positions. Multiply these values by the radius of the outer shell of the head model to determine exact electrode locations.

### Model
- **Head Model**: Spherical three-layer model (refer to Maxwell's equations provided on slides 27-29 for the model details).

### Tasks

- **a)** Generate all possible dipole locations with a resolution of 1 cm and visualize them in 3D. Calculate and save the lead-field matrix.
- **b)** Plot the electrode locations from part (a) on the same 3D plot. Label each electrode with its identifier.
- **c)** Choose a random dipole location (on the sphere's surface) and set its orientation to be radial. Add this dipole’s normalized position and orientation to the plot from part (b).
- **d)** Assign an interictal spike activity from one row of `mat.Interictal` to the chosen dipole. Calculate the potential at 21 electrodes using the dipole’s direction and the lead-field matrix. Plot the electrode potentials over time and label each electrode.
- **e)** Identify the time of the positive peak for all spikes at the electrodes. Define a 7-point window around each peak and store the averaged potentials for all electrodes in a vector. Normalize the potential range (from the most negative to the most positive value) and display this on a color map, marking each electrode's location with a color corresponding to its potential.
- **f)** Repeat part (e) using the function `m.D3_Potential_Display`.

## 2. Inverse Problem Solving

### Tasks

- **g)** Apply MNE and WMNE algorithms to the electrode potentials to solve the inverse problem. 
- **h)** For each localization method used in part (g), estimate the dipole location and orientation. Determine the dipole with the maximum estimated moment and obtain its orientation.
- **i)** Compute the localization error for each method by comparing the estimated location and orientation to the true values.
- **j)** Repeat parts (c) through (i) using a deep dipole model and compare the results with those obtained for the superficial dipole.

## 3. Non-parametric Methods

### Tasks

- **k)** Select two non-parametric methods (e.g., LORETA, LAURA, sLORETA) and repeat parts (g) through (i) using these methods.

## 4. Parametric Model

### Tasks

- **l)** Using a parametric model, solve the inverse problem for a dipole. Define a search space for the dipole location and orientation (e.g., three variables for x, y, z coordinates and/or one variable for the dipole magnitude with a fixed orientation). Use a search algorithm (e.g., genetic algorithm or simulated annealing) to estimate the dipole location and orientation. Repeat parts (g) through (i) for this approach.
- **m)** Choose 15 to 20 adjacent dipoles and set their orientations radially. Plot these dipoles and their orientations as done in part (c).
- **n)** Assign spike activities from `mat.Interictal` to each dipole from part (m). Repeat parts (d), (f), and (g) for these dipoles.
- **o)** For each localization method, compute the estimated dipole moment for each dipole.
- **p)** Plot ROC curves for each method using the results from part (o). Compare the curves, analyzing the overlap between estimated and actual dipoles. Calculate the percentage of correctly and incorrectly identified dipoles.

