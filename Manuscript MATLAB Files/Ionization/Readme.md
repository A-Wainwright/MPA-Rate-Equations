# MATLAB Scripts for Ionization Modeling

## Protein Rate Equation (SRE) Model

This folder contains MATLAB scripts used to simulate strong-field ionization dynamics in protein systems under intense ultrafast laser irradiation.

The model implements a single-rate equation (SRE) framework to describe free-electron generation through strong-field ionization and related nonlinear processes. These simulations were used to generate data and figures presented in the accompanying manuscript and are fully reproducible under the specified conditions.

---

## Contents

This folder contains one simulation script:

* **`ProteinSREModel.m`** — Single-rate equation (SRE) model for strong-field ionization in protein systems

The script is self-contained and includes:

* Laser parameter definitions
* Material parameters
* Strong-field ionization rate calculations
* Time-dependent free-electron density evolution
* Plot generation

---

## Model Description

The simulation is based on a single-rate equation (SRE) formalism describing the evolution of free-electron density under high-intensity femtosecond excitation.

The model includes:

* Strong-field ionization (e.g., Keldysh-type formalism where applicable)
* Avalanche / impact ionization (if enabled)
* Intensity-dependent ionization rates
* Temporal pulse profile modeling

The SRE framework provides a reduced but computationally efficient description of ionization dynamics appropriate for estimating electron densities and assessing damage thresholds in protein systems.

---

## Features

* Simulate strong-field ionization under femtosecond excitation
* Model time-dependent free-electron density
* Estimate ionization thresholds
* Generate manuscript-ready plots
* Fully reproducible simulations

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

