function runDissociation(params, pb)
% pb is optional — all updates are guarded so the function works without it
if nargin < 2, pb = []; end

try
    updatePB(pb, 0.05, 'Validating inputs...');

    % Input validation
    if nargin < 1 || isempty(params)
        error('Parameters structure is required');
    end

    const = getPhysicalConstants();
    % Cross sections
    sigma01 = params.sigma;             % m² (1PA S0→S1)
    sigma02 = params.sigma02;           % m⁴·s/photon (2PA S0→Sn)
    sigma12 = params.sigma12;           % m² (1PA S1→Sn)

    % Time constants (converted to seconds)
    tau_p = params.tau_p;       % Pulse duration: s
    tau_1 = params.tau_1;       % S1 lifetime: s
    tau_2 = params.tau_2;       % S2 lifetime: s
    eps = 1e-15;

    % Laser params
    I0 = params.I * 1e4;                % Peak intensity: W/cm² → W/m²
    wavelength = params.wavelength * 1e-9; % Wavelength: nm → m

    % Sample parameters
    d_sample = params.d_sample * 1e-6;   % sample depth: μm → m

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

    updatePB(pb, 0.12, 'Building time grid...');

    % Set up spatial discretization
    z_max = max(d_sample, 1e-9);
    z = linspace(0, z_max, Nz);
    z_profile = exp(-alpha * z)';  % absorption profile

    % Define pulse region
    pulse_start = -3 * tau_p;
    pulse_end = 3 * tau_p;
    post_pulse_end = max(10 * tau_p, 5 * max(tau_1, tau_2));
    N_pulse = 100;
    t_pulse = linspace(pulse_start, pulse_end, N_pulse);
    N_post = 300;
    post_start = pulse_end + 1e-17;

    % Logarithmic spacing post-pulse
    t_post = logspace(log10(post_start), log10(post_pulse_end), N_post);

    % Full time array (remove duplicates!)
    t = unique([t_pulse, t_post], 'stable');
    Nt = length(t);

    % Initial populations
    n0 = ones(Nz, 1);
    n1 = zeros(Nz, 1);
    n2 = zeros(Nz, 1);
    n0_1 = zeros(Nz, 1);
    n0_2 = zeros(Nz, 1);

    % Precompute the laser pulse profile
    pulse_profile = exp(-4*log(2) * (t/tau_p).^2);
    I_t_array = I0 * pulse_profile;

    % Diagnostics
    max_R01 = I0 * sigma01 / photon_energy;
    max_R02 = (I0^2 / 2) * sigma02 / photon_energy^2;
    max_R12 = I0 * sigma12 / photon_energy;

    fprintf('\nParameters:\n');
    fprintf('  Max R01 rate: %.2e s^-1\n', max_R01);
    fprintf('  Max R02 rate: %.2e s^-1\n', max_R02);
    fprintf('  Max R12 rate: %.2e s^-1\n', max_R12);
    fprintf('  S1 decay rate: %.2e s^-1\n', 1/tau_1);
    fprintf('  S2 decay rate: %.2e s^-1\n', 1/tau_2);

    tic; % timer
    pulse_end_idx = length(t_pulse);

    % ==========================
    % Solver selection
    % ==========================

    % Initial populations
    y0 = [mean(n0); mean(n1); mean(n2); mean(n0_1); mean(n0_2)];
    I_func = @(tt) interp1(t, I_t_array, tt, 'linear', 0)/photon_energy;

    updatePB(pb, 0.20, sprintf('Running %s solver...', params.solver));

    switch params.solver
        case 'Euler'
            fprintf('\nRunning manual Euler solver...\n');

            % Preallocate
            n0_t = zeros(Nt,1); n1_t = zeros(Nt,1); n2_t = zeros(Nt,1);
            n0_1_t = zeros(Nt,1); n0_2_t = zeros(Nt,1);

            % Precompute timesteps
            dt_array = diff([t, t(end)]);

            % Euler loop — update pb every ~10%
            pb_euler_start = 0.20;
            pb_euler_end   = 0.80;

            for i = 1:Nt
                I_t = I_t_array(i);
                dt = dt_array(min(i,end));

                % Rates
                R01 = I_t * sigma01 / photon_energy;
                R12 = I_t * sigma12 / photon_energy;
                R02 = (I_t^2 / 2) * sigma02 / photon_energy^2;

                % Rate equations
                available_01 = max(n0 - n1, 0);
                available_02 = max(n0 - n2, 0);
                available_12 = max(n1 - n2, 0);

                dn0_dt   = -R01*available_01.*z_profile - R02*available_02.*z_profile;
                dn1_dt   = R01*available_01.*z_profile - R12*available_12.*z_profile - n1/tau_1;
                dn2_dt   = R02*available_02.*z_profile + R12*available_12.*z_profile - n2/tau_2;
                dn0_1_dt = n1/tau_1;
                dn0_2_dt = n2/tau_2;

                % Update
                n0   = max(n0   + dn0_dt*dt, 1e-10);
                n1   = max(n1   + dn1_dt*dt, 1e-10);
                n2   = max(n2   + dn2_dt*dt, 1e-10);
                n0_1 = max(n0_1 + dn0_1_dt*dt, 1e-10);
                n0_2 = max(n0_2 + dn0_2_dt*dt, 1e-10);

                % Normalize
                total_pop = n0+n1+n2+n0_1+n0_2;
                mean_total = mean(total_pop);
                if mean_total > 1.0001
                    scale_factor = 1.0/mean_total;
                    n0=n0*scale_factor; n1=n1*scale_factor;
                    n2=n2*scale_factor; n0_1=n0_1*scale_factor; n0_2=n0_2*scale_factor;
                end

                % Store
                n0_t(i)=mean(n0); n1_t(i)=mean(n1); n2_t(i)=mean(n2);
                n0_1_t(i)=mean(n0_1); n0_2_t(i)=mean(n0_2);

                % Update progress bar every ~10% of time steps
                if mod(i, max(1, floor(Nt/10))) == 0 || i == Nt
                    frac = pb_euler_start + (pb_euler_end - pb_euler_start) * i / Nt;
                    phase = 'pulse';
                    if i > pulse_end_idx, phase = 'decay'; end
                    updatePB(pb, frac, sprintf('Euler: %s phase  %.0f%%', phase, 100*i/Nt));
                end

                % Console progress
                if i == pulse_end_idx
                    fprintf('Pulse phase completed (%.1f sec elapsed)\n', toc);
                elseif i == pulse_end_idx + floor((Nt-pulse_end_idx)/2)
                    fprintf('Halfway through decay phase (%.1f sec elapsed)\n', toc);
                elseif mod(i, floor(Nt/5)) == 0
                    fprintf('Progress: %.0f%% (%.1f sec elapsed)\n', 100*i/Nt, toc);
                end
            end

        case 'ode45'
            fprintf('\nRunning ode45 solver...\n');
            odefun = @(tt,n) dissociationODE(tt, n, I_func, sigma01, sigma02, sigma12, tau_1, tau_2, eps);
            [t_sol, y_sol] = ode45(odefun, [t(1) t(end)], y0);
            updatePB(pb, 0.80, 'ode45 complete, unpacking results...');
            t = t_sol;
            n0_t   = y_sol(:,1); n1_t   = y_sol(:,2); n2_t   = y_sol(:,3);
            n0_1_t = y_sol(:,4); n0_2_t = y_sol(:,5);

        case 'ode23'
            fprintf('\nRunning ode23 solver...\n');
            odefun = @(tt,n) dissociationODE(tt, n, I_func, sigma01, sigma02, sigma12, tau_1, tau_2, eps);
            [t_sol, y_sol] = ode23(odefun, [t(1) t(end)], y0);
            updatePB(pb, 0.80, 'ode23 complete, unpacking results...');
            t = t_sol;
            n0_t   = y_sol(:,1); n1_t   = y_sol(:,2); n2_t   = y_sol(:,3);
            n0_1_t = y_sol(:,4); n0_2_t = y_sol(:,5);

        case 'ode23s'
            fprintf('\nRunning ode23s solver...\n');
            odefun = @(tt,n) dissociationODE(tt, n, I_func, sigma01, sigma02, sigma12, tau_1, tau_2, eps);
            [t_sol, y_sol] = ode23s(odefun, [t(1) t(end)], y0);
            updatePB(pb, 0.80, 'ode23s complete, unpacking results...');
            t = t_sol;
            n0_t   = y_sol(:,1); n1_t   = y_sol(:,2); n2_t   = y_sol(:,3);
            n0_1_t = y_sol(:,4); n0_2_t = y_sol(:,5);

        case 'ode15s'
            fprintf('\nRunning ode15s solver...\n');
            odefun = @(tt,n) dissociationODE(tt, n, I_func, sigma01, sigma02, sigma12, tau_1, tau_2, eps);
            [t_sol, y_sol] = ode15s(odefun, [t(1) t(end)], y0);
            updatePB(pb, 0.80, 'ode15s complete, unpacking results...');
            t = t_sol;
            n0_t   = y_sol(:,1); n1_t   = y_sol(:,2); n2_t   = y_sol(:,3);
            n0_1_t = y_sol(:,4); n0_2_t = y_sol(:,5);

        otherwise
            error('Unknown solver: %s', params.solver);
    end

    elapsed_time = toc;
    fprintf('Simulation completed in %.2f seconds\n', elapsed_time);

    updatePB(pb, 0.85, 'Computing final statistics...');

    % Final results
    n_bound_final = n0_t(end) + n1_t(end) + n2_t(end) + n0_1_t(end) + n0_2_t(end);
    n_final_1PA = n0_1_t(end);
    n_final_2PA = n0_2_t(end);
    n_final_total = n_final_1PA + n_final_2PA;

    if n_final_2PA > 0 && n_final_1PA > 0
        pathway_ratio_2PA_to_1PA = n_final_2PA / n_final_1PA;
    else
        pathway_ratio_2PA_to_1PA = 0;
    end

    fprintf('\n=== RESULTS ===\n');
    fprintf('Final populations:\n');
    fprintf('  n0: %.6f, n1: %.6f, n2: %.6f\n', n0_t(end), n1_t(end), n2_t(end));
    fprintf('  n0_1: %.6f, n0_2: %.6f\n', n0_1_t(end), n0_2_t(end));
    fprintf('  Total bound: %.6f\n', n_bound_final);
    fprintf('  2PA/1PA ratio: %.3f\n', pathway_ratio_2PA_to_1PA);

    updatePB(pb, 0.90, 'Generating plots...');

    % Plotting
    plot_data.t = t;
    plot_data.pulse_end_time = pulse_end;
    plot_data.n0_t = n0_t;
    plot_data.n1_t = n1_t;
    plot_data.n2_t = n2_t;
    plot_data.n0_1_t = n0_1_t;
    plot_data.n0_2_t = n0_2_t;

    sampleName = params.sampleName;
    plotDissociationResults(plot_data, params, sampleName);
    updatePB(pb, 0.95, 'Generating pathway comparison plot...');
    plotPathwayComparison(plot_data, params, sampleName);

    updatePB(pb, 1.0, 'Done!');

catch ME
    fprintf('ERROR in runDissociation: %s\n', ME.message);
    fprintf('Error occurred on line %d\n', ME.stack(1).line);
    rethrow(ME);
end
end

% sfetly update 
function updatePB(pb, val, msg)
    if ~isempty(pb) && isvalid(pb)
        pb.Value   = max(0, min(1, val));
        pb.Message = msg;
        drawnow;
    end
end