
function runExcitation(params)
% Timer for the entire function
master_tic = tic;
try
    % Setup / parameter calc
    fprintf('\n=== EXCITATION SIMULATION PARAMETERS ===\n');
    const = getPhysicalConstants();
    wavelength = max(params.wavelength * 1e-9, 1e-10);
    FWHM = max(params.tau_p, const.eps);
    I0 = max(params.I * 1e4, const.eps);
    Mprot = max(params.conc, const.eps);
    molar_extinction_Protein = max(params.alpha, const.eps);
    Cross_section_0_4_per_photon = max(params.sigma02, 0);
    Cross_section_1_4_per_photon = max(params.sigma12, 0);
    t_dephase_1 = max(params.t_dephase_1, const.eps);
    lifetime_1 = max(params.lifetime_1, const.eps);
    lifetime_4 = max(params.lifetime_4, const.eps);
    z_max = max(params.z_max, 1e-9);
    w = 2 * pi * const.c / wavelength;
    E_photon = const.h_bar * w;
    n_prot = Mprot * const.NA * 1e3;
    extinction_coeff = molar_extinction_Protein * Mprot;
    alpha_protein = max(extinction_coeff * 100, const.eps);
    Cross_section_0_1_per_photon = alpha_protein / max(n_prot, const.eps);

    hw = max(const.h_bar * w, const.eps);
    hw2 = max((const.h_bar * w)^2, const.eps);

    sigma_01 = Cross_section_0_1_per_photon / hw;
    sigma_04 = Cross_section_0_4_per_photon / hw2;
    sigma_14 = Cross_section_1_4_per_photon / hw;
    TPA_prot = Cross_section_0_4_per_photon / hw * n_prot;
    t_dephase_4 = 0.7 * t_dephase_1;

    % Parameter Diagnostics
    fprintf('Sample Name: %s\n', params.sampleName);
    fprintf('Wavelength: %.0f nm, Peak Power: %.1e W/cm^2, Pulse Duration: %.0f fs\n', params.wavelength, params.I, FWHM*1e15);
    fprintf('\nParameter Diagnostics:\n');
    fprintf('  Max 1PA Rate (S0->S1): %.2e s^-1\n', sigma_01 * I0);
    fprintf('  Max 2PA Rate (S0->Sn): %.2e s^-1\n', sigma_04 * I0^2);
    fprintf('  Max ESA Rate (S1->Sn): %.2e s^-1\n', sigma_14 * I0);
    fprintf('  S1 Lifetime (τ₁): %.2e s (Rate: %.2e s^-1)\n', lifetime_1, 1/lifetime_1);
    fprintf('  Sn Lifetime (τ₄): %.2e s (Rate: %.2e s^-1)\n', lifetime_4, 1/lifetime_4);
    fprintf('  S1 Dephasing (T₂): %.2e s (Rate: %.2e s^-1)\n', t_dephase_1, 1/t_dephase_1);
    fprintf('  Sn Dephasing (T₂): %.2e s (Rate: %.2e s^-1)\n', t_dephase_4, 1/t_dephase_4);

    % Set up the simulation grid
    z_min = max(0.1e-6, z_max/1000);
    n_depth = min(max(50, 300), 1000);
    z = linspace(z_min, z_max, n_depth);

    t_max = max(2e-12, 10*FWHM);
    t_pre = -2 * FWHM;
    n_time = min(max(50, 100), 500);
    t = linspace(t_pre, t_max, n_time);

    I_t = I0 * exp(-2 * (2*log(2)) * (t / FWHM).^2);
    odeIntensity = @(z_val, I) -max(alpha_protein, const.eps) * max(I, 0) - max(TPA_prot, 0) * max(I, 0).^2;
    try
        [~, I_sol] = ode45(odeIntensity, z, I0, odeset('RelTol', 1e-6, 'AbsTol', 1e-12));
        I_z = max(I_sol', 0);
    catch
        I_z = I0 * exp(-alpha_protein * z);
    end
    I_z_t = zeros(n_depth, n_time);
    for ti = 1:n_time
        I_z_t(:, ti) = I_z * (I_t(ti) / max(I0, const.eps));
    end

    % MAIN SIM LOOP
    fprintf('\nStarting simulation across %d depth points...\n', n_depth);
    sim_tic = tic;

    n0_matrix = zeros(n_depth, n_time); n1_matrix = zeros(n_depth, n_time);
    n1_d_matrix = zeros(n_depth, n_time); n4_matrix = zeros(n_depth, n_time);
    n4_d_matrix = zeros(n_depth, n_time); n01_matrix = zeros(n_depth, n_time);
    n04_matrix = zeros(n_depth, n_time);

    pulse_end_idx = find(t > FWHM*3, 1);
    if isempty(pulse_end_idx), pulse_end_idx = n_time; end

    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, 'NonNegative', 1:7);
    for zi = 1:n_depth
        try
            I_func = @(t_in) interp1(t, I_z_t(zi,:), t_in, 'linear', 0);
            init = [n_prot; 0; 0; 0; 0; 0; 0];

            switch params.solver
                case 'Euler'
                    % Simple fixed-step Euler integrator
                    dt = (t(pulse_end_idx) - t(1)) / (length(t(1:pulse_end_idx)) - 1);
                    n_sol1 = zeros(length(t(1:pulse_end_idx)), length(init));
                    n_sol1(1,:) = init';
                    for k = 2:length(t(1:pulse_end_idx))
                        tn = t(k-1);
                        n_prev = n_sol1(k-1,:)';
                        dn = rate_eqs_during_pulse(tn, n_prev, I_func, ...
                            sigma_01, sigma_04, sigma_14, ...
                            lifetime_1, lifetime_4, t_dephase_1, t_dephase_4, const.eps);
                        n_sol1(k,:) = (n_prev + dt*dn)'; % Euler step
                    end
                    t_sol1 = t(1:pulse_end_idx);

                case 'ode45'
                    ode1 = @(t_in,n) rate_eqs_during_pulse(t_in,n,I_func, ...
                        sigma_01,sigma_04,sigma_14,lifetime_1,lifetime_4, ...
                        t_dephase_1,t_dephase_4,const.eps);
                    [t_sol1, n_sol1] = ode45(ode1, t(1:pulse_end_idx), init, options);

                case 'ode23'
                    ode1 = @(t_in,n) rate_eqs_during_pulse(t_in,n,I_func, ...
                        sigma_01,sigma_04,sigma_14,lifetime_1,lifetime_4, ...
                        t_dephase_1,t_dephase_4,const.eps);
                    [t_sol1, n_sol1] = ode23(ode1, t(1:pulse_end_idx), init, options);

                case 'ode23s'
                    ode1 = @(t_in,n) rate_eqs_during_pulse(t_in,n,I_func, ...
                        sigma_01,sigma_04,sigma_14,lifetime_1,lifetime_4, ...
                        t_dephase_1,t_dephase_4,const.eps);
                    [t_sol1, n_sol1] = ode23s(ode1, t(1:pulse_end_idx), init, options);

                case 'ode15s'
                    ode1 = @(t_in,n) rate_eqs_during_pulse(t_in,n,I_func, ...
                        sigma_01,sigma_04,sigma_14,lifetime_1,lifetime_4, ...
                        t_dephase_1,t_dephase_4,const.eps);
                    [t_sol1, n_sol1] = ode15s(ode1, t(1:pulse_end_idx), init, options);

                otherwise
                    error('Unknown solver: %s', params.solver);
            end

            % After the pulse (same switch structure)
            final_state = n_sol1(end,:)';
            switch params.solver
                case 'Euler'
                    dt = (t(end) - t(pulse_end_idx)) / (length(t(pulse_end_idx:end)) - 1);
                    n_sol2 = zeros(length(t(pulse_end_idx:end)), length(init));
                    n_sol2(1,:) = final_state';
                    for k = 2:length(t(pulse_end_idx:end))
                        tn = t(pulse_end_idx + k - 2);
                        n_prev = n_sol2(k-1,:)';
                        dn = rate_eqs_after_pulse(tn, n_prev, ...
                            lifetime_1,lifetime_4,t_dephase_1,t_dephase_4,const.eps);
                        n_sol2(k,:) = (n_prev + dt*dn)'; % Euler step
                    end
                    t_sol2 = t(pulse_end_idx:end);

                case 'ode45'
                    ode2 = @(t_in,n) rate_eqs_after_pulse(t_in,n, ...
                        lifetime_1,lifetime_4,t_dephase_1,t_dephase_4,const.eps);
                    [t_sol2, n_sol2] = ode45(ode2, t(pulse_end_idx:end), final_state, options);

                case 'ode23'
                    ode2 = @(t_in,n) rate_eqs_after_pulse(t_in,n, ...
                        lifetime_1,lifetime_4,t_dephase_1,t_dephase_4,const.eps);
                    [t_sol2, n_sol2] = ode23(ode2, t(pulse_end_idx:end), final_state, options);

                case 'ode23s'
                    ode2 = @(t_in,n) rate_eqs_after_pulse(t_in,n, ...
                        lifetime_1,lifetime_4,t_dephase_1,t_dephase_4,const.eps);
                    [t_sol2, n_sol2] = ode23s(ode2, t(pulse_end_idx:end), final_state, options);

                case 'ode15s'
                    ode2 = @(t_in,n) rate_eqs_after_pulse(t_in,n, ...
                        lifetime_1,lifetime_4,t_dephase_1,t_dephase_4,const.eps);
                    [t_sol2, n_sol2] = ode15s(ode2, t(pulse_end_idx:end), final_state, options);
            end

            % Concatenate pulse and post-pulse parts
            if exist('t_sol2','var')
                t_sol = [t_sol1; t_sol2(2:end,:)];
                n_sol = [n_sol1; n_sol2(2:end,:)];
            else
                t_sol = t_sol1;
                n_sol = n_sol1;
            end

            % Interpolate results back onto the original time grid
            n0_matrix(zi,:) = interp1(t_sol, n_sol(:,1), t, 'linear', n_prot);
            n1_matrix(zi,:) = interp1(t_sol, n_sol(:,2), t, 'linear', 0);
            n1_d_matrix(zi,:) = interp1(t_sol, n_sol(:,3), t, 'linear', 0);
            n4_matrix(zi,:) = interp1(t_sol, n_sol(:,4), t, 'linear', 0);
            n4_d_matrix(zi,:) = interp1(t_sol, n_sol(:,5), t, 'linear', 0);
            n01_matrix(zi,:) = interp1(t_sol, n_sol(:,6), t, 'linear', 0);
            n04_matrix(zi,:) = interp1(t_sol, n_sol(:,7), t, 'linear', 0);

        catch ME
            warning('Failed at depth index %d: %s', zi, ME.message);
            n0_matrix(zi,:) = n_prot * ones(1, n_time);
        end

        if mod(zi, floor(n_depth/4)) == 0 && zi > 1
            fprintf('  Progress: %.0f%% (%.1f sec elapsed)\n', 100*zi/n_depth, toc(sim_tic));
        end
    end
    fprintf('Simulation completed in %.2f seconds.\n', toc(sim_tic));

    % Post processing and results
    ratio_matrix = n04_matrix ./ max(n01_matrix + const.eps, const.eps);
    z_range = max(z(end) - z(1), const.eps);

    n_prot_safe = max(n_prot, const.eps);
    weighted_avg_1PA_vs_time = trapz(z, n01_matrix / n_prot_safe) / z_range;
    weighted_avg_2PA_vs_time = trapz(z, n04_matrix / n_prot_safe) / z_range;
    weighted_avg_n_0_vs_time = trapz(z, n0_matrix / n_prot_safe) / z_range;
    weighted_avg_n_1_vs_time = trapz(z, n1_matrix / n_prot_safe) / z_range;
    weighted_avg_n_1_d_vs_time = trapz(z, n1_d_matrix / n_prot_safe) / z_range;
    weighted_avg_n_4_vs_time = trapz(z, n4_matrix / n_prot_safe) / z_range;
    weighted_avg_n_4_d_vs_time = trapz(z, n4_d_matrix / n_prot_safe) / z_range;

    final_pop_S0 = weighted_avg_n_0_vs_time(end);
    final_pop_1PA = weighted_avg_1PA_vs_time(end);
    final_pop_2PA = weighted_avg_2PA_vs_time(end);
    total_final_pop = final_pop_S0 + final_pop_1PA + final_pop_2PA;

    fprintf('\n=== SIMULATION RESULTS ===\n');
    fprintf('Final Depth-Averaged Population (%% of initial):\n');
    fprintf('  Remaining Ground State (S₀): %.2f%%\n', final_pop_S0 * 100);
    fprintf('  Relaxed to S₁ State (1PA Product): %.2f%%\n', final_pop_1PA * 100);
    fprintf('  Relaxed to S₄ State (2PA Product): %.2f%%\n', final_pop_2PA * 100);
    fprintf('  Total Accounted For: %.2f%%\n', total_final_pop * 100);

    fprintf('\n=== PATHWAY COMPARISON RESULTS ===\n');
    fprintf('  1-Photon Pathway Yield: %.4f\n', final_pop_1PA);
    fprintf('  2-Photon Pathway Yield: %.4f\n', final_pop_2PA);

    if final_pop_1PA > const.eps
        fprintf('  2PA/1PA Ratio: %.3f\n', final_pop_2PA / final_pop_1PA);
        fprintf('  1PA Contribution: %.1f%%\n', (final_pop_1PA / (final_pop_1PA + final_pop_2PA)) * 100);
        fprintf('  2PA Contribution: %.1f%%\n', (final_pop_2PA / (final_pop_1PA + final_pop_2PA)) * 100);
    else
        fprintf('  1-Photon pathway resulted in negligible yield.\n');
    end

    % Plotting
    fprintf('\nGenerating plots...\n');
    plot_tic = tic;
    try
        fig_pos = [50, 50, 950, 750];
        title_fontsize = 14;
        label_fontsize = 12;
        fluence = I0 * FWHM * 1e-4; % J/cm^2, as in original code

        % Figure 1
        figure('Name', [params.sampleName, ' - Yields & Ratio Map'], 'Position', fig_pos, 'NumberTitle', 'off');
        subplot(3,1,1);
        plot(t*1e15, weighted_avg_1PA_vs_time, 'b', 'LineWidth', 2);
        xlabel('Time (fs)', 'FontSize', label_fontsize);
        ylabel('n_{S1} / n_{prot}', 'FontSize', label_fontsize);
        title(sprintf('%s: %.0f GW/cm^2, %.0f fs, %.3f J/cm^2', params.sampleName, I0/1e9, FWHM*1e15, fluence), 'FontSize', title_fontsize, 'FontWeight', 'bold');
        grid on;

        subplot(3,1,2);
        plot(t*1e15, weighted_avg_2PA_vs_time, 'r', 'LineWidth', 2);
        xlabel('Time (fs)', 'FontSize', label_fontsize);
        ylabel('n_{Sn} / n_{prot}', 'FontSize', label_fontsize);
        grid on;

        subplot(3,1,3);
        imagesc(t*1e15, z*1e6, ratio_matrix);
        xlabel('Time (fs)', 'FontSize', label_fontsize);
        ylabel('Depth (μm)', 'FontSize', label_fontsize);
        title('2PA/1PA Ratio', 'FontSize', title_fontsize, 'FontWeight', 'bold');
        cb = colorbar;
        ylabel(cb, 'Ratio', 'FontSize', label_fontsize);
        axis xy;

        % Figure 2 (focuses on energy level dynamics)
        figure('Name',[params.sampleName, ' - Energy Level Dynamics'], 'Position', fig_pos, 'NumberTitle', 'off');
        ax = gobjects(7,1);
        ax(1) = subplot(7,1,1); plot(t*1e15, weighted_avg_n_0_vs_time, 'b'); ylabel('n_{0}/n_{prot}', 'FontSize', label_fontsize); title(sprintf('%s: %.0f GW/cm^2, %.0f fs', params.sampleName, I0/1e9, FWHM*1e15), 'FontSize', title_fontsize, 'FontWeight', 'bold');
        ax(2) = subplot(7,1,2); plot(t*1e15, weighted_avg_n_1_vs_time, 'r'); ylabel('S1 coh', 'FontSize', label_fontsize);
        ax(3) = subplot(7,1,3); plot(t*1e15, weighted_avg_n_1_d_vs_time, 'r'); ylabel('S1 dephased', 'FontSize', label_fontsize);
        ax(4) = subplot(7,1,4); plot(t*1e15, weighted_avg_n_4_vs_time, 'r'); ylabel('Sn coh', 'FontSize', label_fontsize);
        ax(5) = subplot(7,1,5); plot(t*1e15, weighted_avg_n_4_d_vs_time, 'r'); ylabel('Sn dephased', 'FontSize', label_fontsize);
        ax(6) = subplot(7,1,6); plot(t*1e15, weighted_avg_1PA_vs_time, 'r'); ylabel('S1 relaxed', 'FontSize', label_fontsize);
        ax(7) = subplot(7,1,7); plot(t*1e15, weighted_avg_2PA_vs_time, 'r'); ylabel('Sn relaxed', 'FontSize', label_fontsize);
        xlabel(ax(7), 'Time (fs)', 'FontSize', label_fontsize);
        for i = 1:7, grid(ax(i), 'on'); end

        % Figure 3
        figure('Name',[params.sampleName, ' - S0,S1,Sn Dynamics'], 'Position', fig_pos, 'NumberTitle', 'off');
        ax2 = gobjects(5,1);
        ax2(1) = subplot(5,1,1); plot(t*1e15, weighted_avg_n_0_vs_time, 'b'); ylabel('S0', 'FontSize', label_fontsize); title(sprintf('%s: %.0f GW/cm^2, %.0f fs', params.sampleName, I0/1e9, FWHM*1e15), 'FontSize', title_fontsize, 'FontWeight', 'bold');
        ax2(2) = subplot(5,1,2); plot(t*1e15, weighted_avg_n_1_vs_time + weighted_avg_n_1_d_vs_time, 'r'); ylabel('S1', 'FontSize', label_fontsize);
        ax2(3) = subplot(5,1,3); plot(t*1e15, weighted_avg_n_4_vs_time + weighted_avg_n_4_d_vs_time, 'r'); ylabel('Sn', 'FontSize', label_fontsize);
        ax2(4) = subplot(5,1,4); plot(t*1e15, weighted_avg_1PA_vs_time, 'r'); ylabel('S1 final', 'FontSize', label_fontsize);
        ax2(5) = subplot(5,1,5); plot(t*1e15, weighted_avg_2PA_vs_time, 'r'); ylabel('Sn final', 'FontSize', label_fontsize);
        xlabel(ax2(5), 'Time (fs)', 'FontSize', label_fontsize);
        for i = 1:5, grid(ax2(i), 'on'); end

    catch ME
        warning('Plotting error: %s', ME.message);
    end
    fprintf('Plots generated in %.2f seconds.\n', toc(plot_tic));

catch ME
    rethrow(ME);
end
fprintf('Total runtime for Excitation: %.2f seconds.\n\n', toc(master_tic));
end