function [diff_n_sfi, t_rho_init, delta_tilda_ev_gap] = delta_tilda(gamma_gap, gamma_int, E_gap_eV, E_int_eV, m, w, e, F, t)

    % Function solving for the change in the band gap (Delta tilda) for a given wavelength
    const = getPhysicalConstants();
    
    % Convert energies from eV to Joules
    E_gap = E_gap_eV * const.e; % Convert eV to J
    E_int = E_int_eV * const.e; % Convert eV to J
    
    % Band gap energy modification - vectorized
    MPI_matrix = (gamma_gap > 1);
    TI_matrix = (gamma_gap <= 1);
    
    delta_tilda_TI_gap = 2/pi .* E_gap .* (sqrt(1 + gamma_gap.^2)) ./ gamma_gap .* ...
                         ellipticE(1 ./ (sqrt(1 + gamma_gap.^2))); % full keldysh model in J
    delta_tilda_mpi_gap = E_gap + (const.e^2 .* F.^2) ./ (4 .* const.m .* w.^2); % for the case of MPI from keldysh paper
    delta_tilda_gap = delta_tilda_TI_gap .* TI_matrix + delta_tilda_mpi_gap .* MPI_matrix; % [J]
    delta_tilda_ev_gap = delta_tilda_gap ./ const.e; % [eV]
    
    % Intermediate level energy modification
    MPI_matrix = (gamma_int > 1);
    TI_matrix = (gamma_int <= 1);
    delta_tilda_TI_int = 2/pi .* E_int .* (sqrt(1 + gamma_int.^2)) ./ gamma_int .* ...
                         ellipticE(1 ./ (sqrt(1 + gamma_int.^2))); % full keldysh model in J
    delta_tilda_mpi_int = E_int + (const.e^2 .* F.^2) ./ (4 .* const.m .* w.^2); % for the case of MPI from keldysh paper
    delta_tilda_int = delta_tilda_TI_int .* TI_matrix + delta_tilda_mpi_int .* MPI_matrix; % [J]
    
    % Keldysh strong field ionization 
    rho_ini_max = 1e19 * 1e6; % m^-3 from page 8, converted from cm^-3
    dn_sfi_rho_init_c = n_sfi(delta_tilda_int, gamma_int, w);

    ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    [t_rho_init, rho_init] = ode45(@(t_solve, n) dn_sft_solver(t_solve, n, t, dn_sfi_rho_init_c, rho_ini_max), t, 0, ode_options);

    dn_sfi_rho_gap = n_sfi(delta_tilda_gap, gamma_gap, w);
    
    % Interpolate rho_init to match time vector t
    rho_init_interp = interp1(t_rho_init, rho_init, t, 'linear', 0);
    
    diff_n_sfi = dn_sfi_rho_init_c .* (1 - rho_init_interp / rho_ini_max) + dn_sfi_rho_gap;
end