# Modeling Multiphoton Absorption and Ionization in Ultrafast Pump-Probe Experiments

This repository contains MATLAB code and a standalone GUI for modeling multiphoton absorption processes, used to generate the data presented in the accompanying manuscript. The project is designed both for reproducing simulation results and for optimizing experimental conditions in pump-probe experiments.

## Authors
**Alexander A.C. Wainwright**¹, **Syeda N. Mahdia**², **Khaled Madhoun**¹, **Jessica E. Besaw**³, **R.J. Dwayne Miller**¹ ⁴*  
\* Corresponding author, dwayne.miller@utoronto.ca

1. Dept. of Physics, University of Toronto, Toronto, ON, Canada, M5R 2M8  
2. Dept. of Engineering Science, University of Toronto, Toronto, ON, Canada, M5R 2M8  
3. Dept. of Biochemistry, University of Toronto, Toronto, ON, Canada, M5S 1A8  
4. Dept. of Chemistry, University of Toronto, Toronto, ON, Canada, M5S 3H4

## Contents
- **MATLAB Scripts**: Scripts used to simulate multiphoton absorption and generate the figures in the manuscript.  
- **Executable GUI**: A user-friendly GUI to explore parameter space and optimize experimental conditions.  

## Features

- Simulate one-photon and two-photon absorption in protein and water systems  
- Include depth-dependent power loss and intensity attenuation  
- Optimize pump-probe experimental parameters with a graphical interface  
- Reproduce all figures from the manuscript  

## Installation

### MATLAB Scripts

1. Clone this repository:  
   ```bash
   git clone https://github.com/yourusername/multiphoton-absorption.git
   ```
2. Open MATLAB and navigate to the cloned folder.
3. Run the main script `main_simulation.m` to reproduce the manuscript results.

### GUI

- The GUI is provided as a compiled executable. Simply run the file `MAPS.exe` (Windows) or the corresponding executable for your OS.  
- No MATLAB license is required to run the executable.  

## Usage

### MATLAB Scripts
1. Modify simulation parameters at the top of the `main_simulation.m` script.  
2. Run the script to generate outputs, including plots and data files.

### GUI
1. Launch the executable.  
2. Adjust experimental parameters (laser wavelength, intensity, sample properties, etc.) via the interface.  
3. Press "Run Simulation" to compute absorption profiles and visualize results in real time.  

## Dependencies
- MATLAB R2020b or later (for running scripts)  
   - Symbolic math toolbox
   - Control system toolbox
- The GUI is standalone and requires no MATLAB installation  

## Contributing
Contributions are welcome! Please submit issues or pull requests for bug fixes, feature requests, or improvements.

## License
This project is licensed under the GPL-3.0 license. See the [LICENSE](LICENSE) file for details.

## Citation
If you use this code in your research, please cite the associated manuscript:

Alexander A.C. Wainwright, Syeda N. Mahdia, Khaled Madhoun, Jessica E. Besaw, R.J. Dwayne Miller,Modeling Multiphoton Absorption and Ionization in Ultrafast Pump-Probe Experiments, JCP (Under Review), 2026.
