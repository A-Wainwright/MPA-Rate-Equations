# MATLAB Scripts for Multiphoton Absorption and Ionization Modeling

This folder contains MATLAB scripts used to simulate multiphoton absorption and strong-field ionization processes relevant to pump–probe and high-intensity ultrafast laser experiments. These scripts were used to generate the data and figures presented in the accompanying manuscript.

---

## Contents

This folder contains two subdirectories:

* **`MPA/`** — Multiphoton absorption (MPA) simulations
  Includes rate-equation models for nonlinear absorption processes in systems such as proteins and molecular solutions. These models account for one- and two-photon absorption, saturation effects, and depth-dependent intensity attenuation.

* **`Ionization/`** — Strong-field ionization simulations
  Includes single-rate equation (SRE) models describing free-electron generation under intense femtosecond excitation.

* **`Sensitvity/`** — Sinsitivity analysis of the rate-equation models for nonlinear absorption processes in the MPA folder.

Each subfolder contains self-contained MATLAB scripts, required parameters, and figure-generation routines.

---

## Features

* Simulate one-photon and two-photon absorption in protein and molecular systems
* Model strong-field ionization and free-electron density evolution
* Include depth-dependent power loss and intensity attenuation (MPA models)
* Generate manuscript-ready plots and reproducible data
* Fully reproducible simulations corresponding to reported results

---

## Installation

1. Clone this repository:

   ```bash
   git clone https://github.com/yourusername/multiphoton-absorption.git
   ```

2. Open MATLAB and navigate to this folder.

3. Ensure all subfolders are on the MATLAB path.

---

## Usage

1. Navigate to the desired subfolder:

   * `MPA/` for multiphoton absorption simulations
   * `Ionization/` for strong-field ionization simulations

2. Open the relevant MATLAB script in that folder.

3. Modify simulation parameters at the top of the script (e.g., laser intensity, wavelength, pulse duration, sample properties).

4. Run the script to generate simulation outputs, including plots and data files.

Refer to the README file within each subfolder for system-specific details.

---

## Dependencies
- MATLAB R2020b or later (for running scripts)  
   - Symbolic math toolbox
   - Control system toolbox

---

## Support

For questions regarding the MATLAB scripts, contact:

**R.J. Dwayne Miller**
[dwayne.miller@utoronto.ca](mailto:dwayne.miller@utoronto.ca)

**Alexander A.C. Wainwright**
[Alexander.Wainwright@mail.utoronto.ca](mailto:Alexander.Wainwright@mail.utoronto.ca)
