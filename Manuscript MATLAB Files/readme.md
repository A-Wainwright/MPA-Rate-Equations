# MATLAB Scripts for Multiphoton Absorption Modeling

This folder contains all MATLAB scripts used to simulate multiphoton absorption processes for pump-probe experiments. These scripts were used to generate the data and figures presented in the accompanying manuscript.

## Contents

* **main\_simulation.m** — Main script to run simulations and reproduce manuscript figures.
* **Supporting Scripts** — Additional MATLAB scripts/functions used by `main_simulation.m` for calculations, plotting, and data processing.
* **Data Files** — Any required input data for simulations (e.g., absorption coefficients, experimental parameters).

## Features

* Simulate one-photon and two-photon absorption in protein and water systems
* Include depth-dependent power loss and intensity attenuation
* Generate plots and data files for analysis
* Fully reproducible simulations for manuscript results

## Installation

1. Clone this repository:

   ```bash
   git clone https://github.com/yourusername/multiphoton-absorption.git
   ```
2. Open MATLAB and navigate to the `MATLAB_Files` folder.
3. Ensure all supporting scripts and data files are in the same folder or on the MATLAB path.

## Usage

1. Open `main_simulation.m` in MATLAB.
2. Modify simulation parameters at the top of the script (laser intensity, wavelength, sample properties, etc.).
3. Run the script to generate simulation outputs, including plots and data files.
4. Use supporting scripts for additional analyses if needed.

## Dependencies

* MATLAB R2020b or later
* No additional toolboxes required

## Support

For questions regarding the MATLAB scripts, contact the corresponding author:

**R.J. Dwayne Miller** – dwayne.miller@utoronto.ca
or
**Alexander A.C. Wainwright** – Alexander.Wainwright@mail.utoronto.ca


## Citation
If you use these MATLAB scripts in your research, please cite the associated manuscript:

*Alexander A.C. Wainwright, Syeda Mahdia, Khaled Madhoun, Jessica E. Besaw, R.J. Dwayne Miller, “Title,” Journal, Year.*
