# MATLAB Scripts for Multiphoton Absorption Modeling

## I₃⁻, Bacteriorhodopsin (BR), and Myoglobin (Mb)

This folder contains MATLAB scripts used to simulate multiphoton absorption (MPA) and ionization dynamics in:

* **Bacteriorhodopsin (BR)**
* **Myoglobin (Mb)**
* **Triiodide (I₃⁻)**

These scripts implement coupled rate-equation models to describe nonlinear excitation processes under intense femtosecond laser irradiation. The simulations were used to generate data and figures presented in the accompanying manuscript and are fully reproducible under the specified conditions.

---

## Contents

This folder contains three system-specific simulation scripts:

* **`MPA_Calculations_Br.m`** — Multiphoton absorption modleing for Bacteriorhodopsin (BR).
* **`MPA_Calculations_Mb.m`** — Multiphoton absorption modleing for Myoglobin (Mb).
* **`MPA_Tri_I.m`** — Multiphoton absorption modleing for Triiodide (I₃⁻).

Each script is self-contained and includes:
* Molecular parameters
* Laser input parameters
* Rate-equation implementation
* Plot generation

---

## Model Description

The simulations are based on multi-level rate equation systems describing:

* One-photon absorption (1PA)
* Two-photon absorption (2PA)
* Saturation effects
* Strong-field or effective ionization (where applicable)
* Depth-dependent intensity attenuation (for BR and Mb)

### Three-Level Approximation
Protein systems are modeled using a minimal three-level system:
1. Ground state
2. Excited electronic state
3. Ionized (free-electron) state

This is the simplest physically consistent model capable of capturing nonlinear absorption and saturation behavior.

---

## Features
* Simulate nonlinear absorption under femtosecond excitation
* Model intensity-dependent saturation behavior
* Estimate ionization fractions
* Include depth-dependent power loss (for protein systems)
* Generate manuscript-ready plots directly from each script
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

1. Open the desired system file:

   * `MPA_Calculations_Br.m`
   * `MPA_Calculations_Mb.m`
   * `MPA_Tri_I.m`

2. Modify simulation parameters at the top of the script:

   * Laser wavelength
   * Pulse duration
   * Peak intensity
   * Sample parameters

3. Run the script to generate:

   * Population dynamics
   * Ionization fraction
   * Nonlinear absorption contributions
   * Intensity attenuation profiles (if enabled)

---

## Dependencies

* MATLAB R2020b or later
* No additional toolboxes required

---

## Limitations

* Reduced-level electronic structure model
* Cross-sections treated as constant unless specified
* No thermal or structural dynamics included
* Plasma effects and beam reshaping not fully modeled
* Parameters derived from literature values and may vary between sources

Results should be interpreted as mechanistic and order-of-magnitude estimates rather than ab initio predictions.

---
## Support

For questions regarding the MATLAB scripts, contact:

**R.J. Dwayne Miller** – [dwayne.miller@utoronto.ca](mailto:dwayne.miller@utoronto.ca)
or
**Alexander A.C. Wainwright** – [Alexander.Wainwright@mail.utoronto.ca](mailto:Alexander.Wainwright@mail.utoronto.ca)
