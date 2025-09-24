function runDissociation(params)
    try
        % Input validation
        if nargin < 1 || isempty(params)
            error('Parameters structure is required');
        end
        
        const = getPhysicalConstants();
        % Cross sections
        sigma01 = params.sigma01;           % m² (1PA S0→S1)
        sigma02 = params.sigma02;           % m⁴·s/photon (2PA S0→S2)
        sigma12 = params.sigma12;           % m² (1PA S1→S2)
        
        % Time constants (converted to seconds)
        tau_p = params.tau_p * 1e-15;       % Pulse duration: fs → s
        tau_1 = params.tau_1 * 1e-15;       % S1 lifetime: fs → s  
        tau_2 = params.tau_2 * 1e-15;       % S2 lifetime: fs → s
        
        % Laser params
        I0 = params.I * 1e4;                % Peak intensity: W/cm² → W/m²
        wavelength = params.wavelength * 1e-9; % Wavelength: nm → m
        
        % Sample parameters
        d_sample = params.d_sample * 1e-6;   % sample depth: μm → m
        
        if isfield(params, 'k_dissoc_1') && ~isempty(params.k_dissoc_1)
            k_dissoc_1 = params.k_dissoc_1; % Dissociation rate from S1 (s^-1)
        else
            k_dissoc_1 = 1e12;  % fallback value - default dissociation value to S1
        end
        
        if isfield(params, 'k_dissoc_2') && ~isempty(params.k_dissoc_2)
            k_dissoc_2 = params.k_dissoc_2; % Dissociation rate from S2 (s^-1)
        else
            k_dissoc_2 = 5e12;  % set up default dissociation value to S2
        end
        
        if isfield(params, 'alpha') && ~isempty(params.alpha)
            alpha = params.alpha;
        else
            alpha = 1e4;  
        end
        
        if isfield(params, 'Nz') && ~isempty(params.Nz)
            Nz = max(10, min(params.Nz, 100));
        else
            Nz = 20; 
        end

        % Constants
        photon_energy = const.h_bar * 2*pi*const.c / wavelength;
        
        % Set up spatial discretization
        z_max = max(d_sample, 1e-9);
        z = linspace(0, z_max, Nz);
        z_profile = exp(-alpha * z)';
        
        % Define pulse region
        pulse_start = -3 * tau_p;  
        pulse_end = 3 * tau_p;     
        post_pulse_end = max(10 * tau_p, 5 * max(tau_1, tau_2));
        N_pulse = 100;
        t_pulse = linspace(pulse_start, pulse_end, N_pulse);
        N_post = 300;  % Number of post-pulse points to evaluate at 
        post_start = pulse_end + 1e-17;  % avoid log(0), make calculations faster 

        % Geometrically space out post-pulse times
        t_post = logspace(log10(post_start), log10(post_pulse_end), N_post);

        % Ensure the post-pulse time starts right after the pulse end
        t_post = t_post(t_post > pulse_end);
        t_post = [pulse_end, t_post]; 
        
        % Combine full time array
        t = [t_pulse, t_post];
        Nt = length(t);
        
        % Initial populations - all in ground state
        n0 = ones(Nz, 1);
        n1 = zeros(Nz, 1);
        n2 = zeros(Nz, 1);
        n0_1 = zeros(Nz, 1);  % Ground state from S1 decay (non-dissociated)
        n0_2 = zeros(Nz, 1);  % Ground state from S2 decay (non-dissociated)
        dissociated = zeros(Nz, 1); 
        
        % Pre-allocate storage arrays
        n0_t = zeros(Nt, 1);
        n1_t = zeros(Nt, 1);
        n2_t = zeros(Nt, 1);
        n0_1_t = zeros(Nt, 1);
        n0_2_t = zeros(Nt, 1);
        dissociated_t = zeros(Nt, 1);
        
        % Precompute the laser pulse profile 
        pulse_profile = exp(-4*log(2) * (t/tau_p).^2);
        I_t_array = I0 * pulse_profile;
        
        % Precompute timesteps so its faster 
        dt_array = zeros(Nt-1, 1);
        for i = 1:(Nt-1)
            dt_array(i) = t(i+1) - t(i);
        end
        
        % Print to the cmd window for diagnostics 
        max_R01 = I0 * sigma01 / photon_energy;
        max_R02 = (I0^2 / 2) * sigma02 / photon_energy^2;
        max_R12 = I0 * sigma12 / photon_energy;
        
        fprintf('\nParameters:\n');
        fprintf('  Max R01 rate: %.2e s^-1\n', max_R01);
        fprintf('  Max R02 rate: %.2e s^-1\n', max_R02);
        fprintf('  Max R12 rate: %.2e s^-1\n', max_R12);
        fprintf('  S1 decay rate: %.2e s^-1\n', 1/tau_1);
        fprintf('  S2 decay rate: %.2e s^-1\n', 1/tau_2);
        fprintf('  S1 dissociation rate: %.2e s^-1\n', k_dissoc_1);
        fprintf('  S2 dissociation rate: %.2e s^-1\n', k_dissoc_2);
        
        % Main simulation loop 
        tic; % Start timer

        % Progress tracking with different phases
        pulse_end_idx = N_pulse;
        
        for i = 1:Nt
            I_t = I_t_array(i);
            
            if i < Nt
                dt = dt_array(i);
            else
                dt = dt_array(end);  
            end
            
            % Calculate transition rates
            R01 = I_t * sigma01 / photon_energy;
            R12 = I_t * sigma12 / photon_energy;
            R02 = (I_t^2 / 2) * sigma02 / photon_energy^2;
            
            % Population dynamics section 
            % Ensure no populations are negative
            available_01 = max(n0 - n1, 0);
            available_02 = max(n0 - n2, 0);
            available_12 = max(n1 - n2, 0);
            
            % Set up rate equations 
            % Ground state: loses population via 1PA and 2PA
            dn0_dt = -R01 * available_01 .* z_profile - R02 * available_02 .* z_profile;
            
            % S1 state: gains from 1PA, loses to S2 via 1PA should decay
            dn1_dt = R01 * available_01 .* z_profile - R12 * available_12 .* z_profile - n1 / tau_1;
            
            % S2 state: gains from 2PA and 1PA from S1, decays naturally  
            dn2_dt = R02 * available_02 .* z_profile + R12 * available_12 .* z_profile - n2 / tau_2;
            
            % Following your original approach - n0_1 and n0_2 track total decay
            dn0_1_dt = n1 / tau_1;
            dn0_2_dt = n2 / tau_2;
            
            ddissociated_dt = (n1 * k_dissoc_1 + n2 * k_dissoc_2);
            
            % Set up adaptive timestep integration 
            n0 = max(n0 + dn0_dt * dt, 1e-10);
            n1 = max(n1 + dn1_dt * dt, 1e-10);
            n2 = max(n2 + dn2_dt * dt, 1e-10);
            n0_1 = max(n0_1 + dn0_1_dt * dt, 1e-10);
            n0_2 = max(n0_2 + dn0_2_dt * dt, 1e-10);
            dissociated = max(dissociated + ddissociated_dt * dt, 1e-10);
            total_pop = n0 + n1 + n2 + n0_1 + n0_2 + dissociated;
            mean_total = mean(total_pop);
            if mean_total > 1.0001  
                scale_factor = 1.0 / mean_total;
                n0 = n0 * scale_factor;
                n1 = n1 * scale_factor;
                n2 = n2 * scale_factor;
                n0_1 = n0_1 * scale_factor;
                n0_2 = n0_2 * scale_factor;
                dissociated = dissociated * scale_factor;
            end
        
            % Store data 
            n0_t(i) = mean(n0);
            n1_t(i) = mean(n1);
            n2_t(i) = mean(n2);
            n0_1_t(i) = mean(n0_1);
            n0_2_t(i) = mean(n0_2);
            dissociated_t(i) = mean(dissociated);
            
            % Progress indicator - prints to the cmd window 
            if i == pulse_end_idx
                fprintf('Pulse phase completed (%.1f sec elapsed)\n', toc);
            elseif i == pulse_end_idx + floor(length(t_post)/2)
                fprintf('Halfway through decay phase (%.1f sec elapsed)\n', toc);
            elseif mod(i, floor(Nt/5)) == 0
                fprintf('Progress: %.0f%% (%.1f sec elapsed)\n', 100*i/Nt, toc);
            end
        end
        
        elapsed_time = toc;
        fprintf('Simulation completed in %.2f seconds\n', elapsed_time);
        
        % Calculate final results
        n_bound_final = n0_t(end) + n1_t(end) + n2_t(end) + n0_1_t(end) + n0_2_t(end);
        dissociation_yield = dissociated_t(end);
        total_accounted = n_bound_final + dissociation_yield;
        
        % Pathway analysis section
        n_final_1PA = n0_1_t(end);
        n_final_2PA = n0_2_t(end);
        n_final_total = n_final_1PA + n_final_2PA;
        
        if n_final_2PA > 0 && n_final_1PA > 0
            pathway_ratio_2PA_to_1PA = n_final_2PA / n_final_1PA;
        else
            pathway_ratio_2PA_to_1PA = 0;
        end
        
        fprintf('\n=== SIMULATION RESULTS ===\n');
        fprintf('Final populations:\n');
        fprintf('  Ground state (original): %.6f\n', n0_t(end));
        fprintf('  S1 state: %.6f\n', n1_t(end));
        fprintf('  S2 state: %.6f\n', n2_t(end));
        fprintf('  Ground state (from S1): %.6f\n', n0_1_t(end));
        fprintf('  Ground state (from S2): %.6f\n', n0_2_t(end));
        fprintf('  Total bound: %.6f\n', n_bound_final);
        fprintf('  Dissociated: %.6f (%.1f%%)\n', dissociation_yield, dissociation_yield * 100);
        fprintf('  Total accounted: %.6f\n', total_accounted);
        
        % Peak populations during pulse vs after pulse
        pulse_indices = 1:pulse_end_idx;
        post_indices = (pulse_end_idx+1):Nt;
        
        fprintf('\nPeak populations:\n');
        fprintf('  During pulse - Max S1: %.6f, Max S2: %.6f\n', ...
                max(n1_t(pulse_indices)), max(n2_t(pulse_indices)));
        if ~isempty(post_indices)
            fprintf('  After pulse - Max S1: %.6f, Max S2: %.6f\n', ...
                    max(n1_t(post_indices)), max(n2_t(post_indices)));
        end
        
        fprintf('\n=== PATHWAY COMPARISON RESULTS ===\n');
        fprintf('Final populations:\n');
        fprintf('  1-Photon pathway (n_0_1): %.6f\n', n_final_1PA);
        fprintf('  2-Photon pathway (n_0_2): %.6f\n', n_final_2PA);
        fprintf('  Total final population: %.6f\n', n_final_total);
        fprintf('  2PA/1PA ratio: %.3f\n', pathway_ratio_2PA_to_1PA);
        
        if n_final_total > 0
            fprintf('  1PA contribution: %.1f%%\n', (n_final_1PA/n_final_total)*100);
            fprintf('  2PA contribution: %.1f%%\n', (n_final_2PA/n_final_total)*100);
        end
        
        fprintf('\nGenerating plots...\n');
        tic;
        
        % Pass adaptive time grid and pulse boundary to plotting functions
        plot_data.t = t;
        plot_data.pulse_end_time = pulse_end;
        plot_data.n0_t = n0_t;
        plot_data.n1_t = n1_t;
        plot_data.n2_t = n2_t;
        plot_data.n0_1_t = n0_1_t;
        plot_data.n0_2_t = n0_2_t;
        plot_data.dissociated_t = dissociated_t;
        
        sampleName = params.sampleName; % Passed the sample name for plots 
        plotDissociationResults(plot_data, params, sampleName);
        plotPathwayComparison(plot_data, params, sampleName);
        
        fprintf('Plots generated in %.2f seconds\n', toc);
        fprintf('Total runtime: %.2f seconds\n', elapsed_time + toc);
    catch ME
        fprintf('ERROR in runDissociation: %s\n', ME.message);
        fprintf('Error occurred on line %d\n', ME.stack(1).line);
        rethrow(ME);
    end
end