# MATLAB Scripts for Sensitivity Analysis of MPA

This folder contains MATLAB scripts used to model the sensitivity of the multiple rate equation models used.

---

## Contents

This folder contains one simulation script:

* **`Depth_MPA_Calculations_Br_Sensitivity_Final_Heatmap.m'** — Sensitivity model used to generate heat map shown in main manuscript
* **Depth_MPA_Calculations_Mb_Sensitivity_Final_V2.m** — Sensitivity model used to generate the sensitivity analysis in the supporting documentation

The script is self-contained and includes:

* Laser parameter definitions
* Material parameters
* Strong-field ionization rate calculations
* Time-dependent free-electron density evolution
* Plot generation

---

## Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/multiphoton-absorption.git
   ```

2. Open MATLAB and navigate to this folder.

3. Ensure the folder is on the MATLAB path.

---

## Usage

1. Open `ProteinSREModel.m` in MATLAB.

2. Modify simulation parameters at the top of the script:

   * Laser wavelength
   * Pulse duration
   * Peak intensity
   * Material parameters

3. Run the script to generate:

   * Time-dependent electron density
   * Ionization rates
   * Threshold estimates
   * Associated plots

---

## Dependencies

* MATLAB R2020b or later
* No additional toolboxes required

---

## Limitations

* Reduced single-rate description (no explicit multi-level population tracking)
* Material parameters derived from literature estimates
* No spatial beam propagation or plasma-induced beam reshaping
* Thermal and structural dynamics not included

Results should be interpreted as mechanistic and order-of-magnitude estimates rather than full ab initio ionization simulations.

---

## Support

For questions regarding the MATLAB scripts, contact:

**R.J. Dwayne Miller**
[dwayne.miller@utoronto.ca](mailto:dwayne.miller@utoronto.ca)

**Alexander A.C. Wainwright**
[Alexander.Wainwright@mail.utoronto.ca](mailto:Alexander.Wainwright@mail.utoronto.ca)

