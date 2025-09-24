function const = getPhysicalConstants()
    const.c = 2.99792458e8;               % Speed of light (m/s)
    const.e = 1.60218e-19;                % Elementary charge (C)
    const.m = 9.11e-31;                   % Electron mass (kg)
    const.m_c = const.m / 2;              % Effective mass in CB
    const.h_bar = 1.054571817e-34;        % Reduced Planck's constant (J·s)
    const.h_bar_ev= 6.5821e-16;           % Planck's reduced Constant in eV·s
    const.epsilon0 = 8.854e-12;           % Vacuum permittivity (F/m)
    const.epsilon0_ev=55.26349406e6;      % Permittivity of free space in e²·eV⁻¹·m⁻¹
    const.NA = 6.022e23;                  % Avogadro's number (1/mol)
    const.eps = 1e-12;                    % Small epsilon to prevent division by zero
    const.M_h20 = 29.9e-27;               % Molar mass of water   
    const.mu0=4*pi*1e-7;                  % Permeability of free space 
end