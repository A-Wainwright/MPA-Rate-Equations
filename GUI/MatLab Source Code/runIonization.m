
function results = runIonization(params)
    const = getPhysicalConstants();
    % Start timer
    tic;
    try
        % Input validation and defaults
        if nargin < 1 || isempty(params)
            params.wavelength = 515e-9; % m
            params.power = 6e12; % W/cm^2 FWHM
            params.tau_p = 250e-15; % s
            params.n_time = 400;
        end
        
        % Set up default for n_time if not provided
        if ~isfield(params, 'n_time') || isempty(params.n_time)
            params.n_time = 400;
        end
        % Set defaults for missing fields
        if ~isfield(params, 'Peak_power_FWHM')
            params.Peak_power_FWHM = params.power;
        end
        if ~isfield(params, 'E_gap')
            params.E_gap = 6.5; % eV, typical for water
        end
        if ~isfield(params, 'E_int')
            params.E_int = 3.2; % eV, intermediate level
        end
        if ~isfield(params, 't_coll')
            params.t_coll = 1e-15; % s, collision time
        end
        if ~isfield(params, 'n_bound')
            params.n_bound = 3.3e28; % m^-3, bound electron density
        end
        if ~isfield(params, 'sampleName')
            params.sampleName = 'Sample'; % Default sample name
        end
        
        % Print input parameters
        fprintf('\n=== Simulation Parameters===\n');
        fprintf('Wavelength: %.0f nm\n', params.wavelength*1e9);
        fprintf('Peak Power (FWHM): %.2f TW/cm^2\n', params.Peak_power_FWHM/1e12);
        fprintf('Pulse Duration: %.0f fs\n', params.tau_p*1e15);
        fprintf('Band Gap Energy: %.2f eV\n', params.E_gap);
        fprintf('Intermediate Level: %.2f eV\n', params.E_int);
        fprintf('Total Time Points: %d\n', params.n_time);
        fprintf('========================================\n\n');
        
        % Set up effective mass 
        m_prime = 0.5 * const.m; % effective mass of electron 
        m_c = const.m / 2; % mass of electron in conduction band 
        % Set up laser intensities - pre-compute constants
        I_max_cm_1 = 2*sqrt(log(2))/sqrt(pi) * params.Peak_power_FWHM;
        I_max_1 = I_max_cm_1 * 1e4; % incident laser intensity W/m^2
        T_laser_1 = params.tau_p; % s IR
        % OPTIMIZED TIME GRID - vectorized creation
        length_measure = 2 * T_laser_1;
        pulse_width = T_laser_1;
        
        % Create adaptive time grid more efficiently
        t_pulse_region = 2 * pulse_width;
        n_pulse = round(0.7 * params.n_time);
        n_outer = params.n_time - n_pulse;
        
        % Vectorized time grid creation
        t_pulse = linspace(-t_pulse_region, t_pulse_region, n_pulse);
        t_left = linspace(-length_measure, -t_pulse_region, round(n_outer/2));
        t_right = linspace(t_pulse_region, length_measure, round(n_outer/2));
        
        % Combine and sort efficiently
        t = [t_left, t_pulse, t_right];
        t = unique(sort(t)); % Combined sort and unique
              
        % Pre-compute time-dependent intensity (vectorized)
        t_norm = (2 * t / T_laser_1).^2;
        I_1 = I_max_1 * exp(-log(2) * t_norm); % more efficient than 2^(-x)
        % Calculate refractive index for water (vectorized)
        wavelength_um_1 = params.wavelength * 1e6; % um
        wl2 = wavelength_um_1^2; % pre-compute wavelength squared
        
        % Optimized Sellmeier equation for water
        n_0_1_sq = 1 + 5.666959820e-1/(1-5.084151894e-3/wl2) + ...
                      1.731900098e-1/(1-1.818488474e-2/wl2) + ...
                      2.095951857e-2/(1-2.625439472e-2/wl2) + ...
                      1.125228406e-1/(1-1.073842352e1/wl2);
        n_0_1 = sqrt(n_0_1_sq);
        fprintf('Status: Computing rate equation constants...\n');
        
        % Pre-compute common constants
        w_1 = 2 * pi * const.c / params.wavelength; % angular frequency
        sqrt_factor = sqrt(2 / (n_0_1 * const.epsilon0 * const.c));
        F_1 = sqrt_factor * sqrt(I_1); % electric Field Intensity in V/m
        
        % Pre-compute Keldysh parameter components
        keldysh_factor_gap = w_1 / const.e * sqrt(m_prime * const.c * n_0_1 * const.epsilon0 * params.E_gap * const.e);
        keldysh_factor_int = w_1 / const.e * sqrt(m_prime * const.c * n_0_1 * const.epsilon0 * params.E_int * const.e);
        
        gama_gap_1 = keldysh_factor_gap ./ sqrt(I_1); % Keldysh parameter
        gama_int_1 = keldysh_factor_int ./ sqrt(I_1); % Keldysh parameter
        fprintf('Status: Computing delta tilda and SFI rates...\n');
        
        % Perform delta_tilda calculation - now using your helper function
        [diff_n_sfi_1, ~, delta_tilda_ev_gap_1] = delta_tilda(gama_gap_1, gama_int_1, params.E_gap, params.E_int, const.m, w_1, const.e, F_1, t);
        
        % Pre-compute avalanche ionization constants
        delta_tilda_gap_1 = delta_tilda_ev_gap_1 * 1.60218e-19;
        w_tc_factor = 1 / (w_1^2 * params.t_coll^2 + 1);
        ai_constant = log(2) * params.t_coll * w_tc_factor * const.e^2 / (const.c * n_0_1 * const.epsilon0 * m_c * 1.5);
        
        n_AI_Asym_1 = ai_constant * I_1 ./ delta_tilda_gap_1;
        fprintf('Status: Solving rate equations...\n');

        ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', max(diff(t))/2, ...
                           'Vectorized', 'on', 'JPattern', 1);
        % Pre-interpolate to solve the ODE faster 
        t_interp = t;
        diff_n_sfi_interp = diff_n_sfi_1(1,:);
        n_AI_interp = n_AI_Asym_1;
        
        % Total number of electrons (SRE1) - using your optimized helper function
        [t_n_total, n_total] = ode45(@(t_solve,n) SRE1(t_solve, n, t_interp, diff_n_sfi_interp, n_AI_interp, params.n_bound), t, 0, ode_options);
        n_total_cm_3 = n_total * 1e-6; % Vectorized conversion
        
        % SFI electrons only - optimized
        [t_SFI, n_SFI] = ode45(@(t_sln,n) interp1(t, diff_n_sfi_interp, t_sln, 'linear', 0), t, 0, ode_options);
        n_SFI_cm_3 = n_SFI * 1e-6; % Vectorized conversion
        
        % AI electrons (difference) - optimized interpolation
        n_SFI_cm_3_interp = interp1(t_SFI, n_SFI_cm_3, t_n_total, 'linear', 0);
        n_AI_cm_3 = n_total_cm_3 - n_SFI_cm_3_interp;
        fprintf('Status: Computing final statistics...\n');
        
        % Calculate final electron densities and percentages
        final_total = max(n_total_cm_3);
        final_SFI = max(n_SFI_cm_3);
        final_AI = max(n_AI_cm_3);
        
        sfi_percentage = (final_SFI / final_total) * 100;
        ai_percentage = (final_AI / final_total) * 100;
        
        % Determine dominant ionization mechanism
        if sfi_percentage > ai_percentage
            dominant_mechanism = 'Strong Field Ionization (SFI)';
        else
            dominant_mechanism = 'Avalanche Ionization (AI)';
        end
        
        % Print final results
        fprintf('\n=== SIMULATION RESULTS ===\n');
        fprintf('Final Total Electron Density: %.2e cm^-3\n', final_total);
        fprintf('Final SFI Electron Density: %.2e cm^-3 (%.1f%%)\n', final_SFI, sfi_percentage);
        fprintf('Final AI Electron Density: %.2e cm^-3 (%.1f%%)\n', final_AI, ai_percentage);
        fprintf('Dominant Mechanism: %s\n', dominant_mechanism);
        
        % Calculate critical density and other metrics
        n_crit = w_1^2 * m_c * const.epsilon0 / const.e^2 * 1e-6; % cm^-3
        fprintf('Critical Density: %.2e cm^-3\n', n_crit);
        
        if final_total > n_crit
            fprintf('STATUS: Plasma formation threshold EXCEEDED\n');
        else
            fprintf('STATUS: Plasma formation threshold NOT reached\n');
        end
        
        fprintf('==========================\n');
        
        % Store results in structure
        results = struct();
        results.t = t_n_total;
        results.n_total_cm3 = n_total_cm_3;
        results.n_SFI_cm3 = n_SFI_cm_3_interp;
        results.n_AI_cm3 = n_AI_cm_3;
        results.final_total = final_total;
        results.final_SFI = final_SFI;
        results.final_AI = final_AI;
        results.sfi_percentage = sfi_percentage;
        results.ai_percentage = ai_percentage;
        results.dominant_mechanism = dominant_mechanism;
        results.n_crit = n_crit;
        results.params = params;
        
        % Plotting
        fprintf('Status: Generating plots...\n');

        figure('Name', [params.sampleName ' - Ionization Results'], 'NumberTitle', 'off');
        clf;
        
        % Pre-convert time to fs for plotting
        t_fs = t_n_total * 1e15;
        semilogy(t_fs, n_total_cm_3, 'Color', [0.0 0.45 0.74], 'LineWidth', 1.25, 'DisplayName', 'Total Ionization'); 
        hold on;
        semilogy(t_fs, n_AI_cm_3, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.25, 'DisplayName', sprintf('AI (%.1f%%)', ai_percentage)); 
        semilogy(t_fs, n_SFI_cm_3_interp, 'Color', [0.47 0.67 0.19], 'LineWidth', 1.25, 'DisplayName', sprintf('SFI (%.1f%%)', sfi_percentage));
        ylim([10^10, 10^22]);
        xlim([-500, 500]);
        
        title(sprintf('%s: %.2f TW/cm^2, %.0f nm\n%s Dominant', ...
              params.sampleName, params.Peak_power_FWHM/1e12, params.wavelength*1e9, dominant_mechanism));
        
        ylabel('Electron density (cm^{-3})');
        xlabel('Time (fs)');
        legend('Location', 'southeast');
        grid on;
        
        computation_time = toc;
        fprintf('Total computation time: %.2f seconds\n', computation_time);
        results.computation_time = computation_time;
        
    catch ME
        fprintf('ERROR in runIonization: %s\n', ME.message);
        fprintf('Error occurred on line %d\n', ME.stack(1).line);
        rethrow(ME);
    end
end
