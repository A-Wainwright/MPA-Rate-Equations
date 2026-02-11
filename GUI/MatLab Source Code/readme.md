# Multiphoton Absorption GUI

This folder contains the MATLAB source code for a graphical user interface (GUI) designed to explore and optimize multiphoton absorption (MPA), excitation, dissociation, and ionization simulations used in ultrafast pump–probe experiments.

The GUI provides an interactive environment for adjusting experimental parameters and visualizing nonlinear absorption and ionization dynamics.

---

## Contents

### Main GUI File

* **`MRE_GUI_Main.m`** — Main entry point for launching the GUI

### Core Modeling Files

* **`SRE1.m`** — Single-rate equation (SRE) ionization model
* **`dn_sfi_solver.m`** — Strong-field ionization solver
* **`n_sfi.m`** — Strong-field ionization rate calculations
* **`dissociationODE.m`** — ODE system for dissociation modeling
* **`rate_eqs_during_pulse.m`** — Rate equations during pulse excitation
* **`rate_eqs_after_pulse.m`** — Rate equations following pulse excitation

### Simulation and Execution Scripts

* **`runExcitation.m`** — Excitation simulations
* **`runDissociation.m`** — Dissociation simulations
* **`runIonization.m`** — Ionization simulations

### Plotting Utilities

* **`plotDissociationResults.m`**
* **`plotPathwayComparison.m`**

### Supporting Functions

* **`Q.m`**
* **`delta_tilda.m`**
* **`getPhysicalConstants.m`**
* **`checkCompilationReadiness.m`**

---

## Features

* Interactive adjustment of:

  * Laser wavelength
  * Pulse duration
  * Peak intensity
  * Material and molecular parameters

* Simulation of:

  * Multiphoton excitation
  * Strong-field ionization
  * Dissociation pathways
  * Free-electron density evolution

* Real-time visualization of:

  * Population dynamics
  * Ionization fractions
  * Pathway comparisons

* Designed for both interactive exploration and reproducible manuscript simulations

---

## Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/multiphoton-absorption.git
   ```

2. Open MATLAB (R2020b or later recommended).

3. Navigate to this folder.

4. Ensure all files are on the MATLAB path.

---

## Usage

1. Open MATLAB.

2. Run:

   ```matlab
   MRE_GUI_Main
   ```

3. Adjust parameters within the GUI interface.

4. Run simulations and visualize results directly within the GUI.

For compiled standalone versions, ensure all required MATLAB Runtime dependencies are installed.

---

## Requirements

* MATLAB R2020b or later
* No additional toolboxes required (unless compiling a standalone version)

---

## Support

For questions or issues related to the GUI, please contact:

**R.J. Dwayne Miller**
[dwayne.miller@utoronto.ca](mailto:dwayne.miller@utoronto.ca)

**Alexander A.C. Wainwright**
[Alexander.Wainwright@mail.utoronto.ca](mailto:Alexander.Wainwright@mail.utoronto.ca)
