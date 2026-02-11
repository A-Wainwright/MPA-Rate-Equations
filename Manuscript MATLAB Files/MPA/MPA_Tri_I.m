%% MPA modleing of I3- Under strong and normal  
% Alexander A.C. Wainwright
clear; clc; close all;

%% --- Constants and Parameters ---
c = 299792458; % speed of light in m/s
e = 1.60217663e-19; % charge of the electron in C
h_bar = 1.054571817e-34; % reduced Planck's constant in J/Hz

% Wavelength and angular frequency
wavelength = 400e-9; % [m]
w = 2 * pi * c / wavelength;
E_photon = h_bar * w; % photon energy [J]

% Protein density and cross sections
M_prot_solution = 2.9060; % [M]
NA = 6.022e23; % Avogadro's number
n_prot = M_prot_solution * NA * 1e3; % [1/m^3]
molar_extinction_Protein = 5000; % M-1 cm-1 at 400 nm

% Cross sections calculation
extinction_coefficient_proteins = molar_extinction_Protein * M_prot_solution; % [1/cm]
alpha_proteins = extinction_coefficient_proteins *100; % [1/m]
Cross_section_1 = alpha_proteins/n_prot; % [m^2]

% Experimental cross sections
Cross_section_0_1_per_photon = Cross_section_1; % [m^2]
Cross_section_0_4_per_photon = 100e-58; % [m^4*s/photon]
Cross_section_1_4_per_photon = Cross_section_1*0.1; % [m^2]

% Convert to intensity-based cross sections
sigma_01 = Cross_section_0_1_per_photon / (h_bar * w); % [m^2 / (W/m^2)]
sigma_04 = Cross_section_0_4_per_photon / (h_bar * w)^2; % [m^4 / (W/m^2)^2]
sigma_14 = Cross_section_1_4_per_photon / (h_bar * w); % [m^2 / (W/m^2)]

% Two-photon absorptivity [1/(m·W)]
TPA_prot = Cross_section_0_4_per_photon / (h_bar * w) * n_prot;

% New time constants from the provided values
lifetime_1 = 100e-15; %[s]

lifetime_4 =0.7*lifetime_1; %[s]

% Laser pulse parameters
FWHM = 100e-15; %[s]
Peak_power_FWHM = 1.67e9; %[W/cm^2]
I0_cm = Peak_power_FWHM; %[W/cm^2]
I0 = I0_cm * 1e4; %[W/m^2]

% Time and depth grid
t_max = 0.5e-12; %Maximum time [s]
t_pre = -2 * FWHM; %Start time [s]
n_time = 100; %Number of time points
z_max = 120e-9; %Maximum depth [m]
z_min = 0.1e-6; %Minimum depth [m]
n_depth = 300; %Number of depth points

t = linspace(t_pre, t_max, n_time);
z = linspace(z_min, z_max, n_depth);

% Pre-calculate the time-dependent intensity at the surface
I_t = I0 * 2.^(-(2 * (t / FWHM)).^2);

%% --- Solve for Depth-Dependent Intensity using ODE45 ---
% Define ODE for intensity attenuation (one-photon + two-photon absorption)
odeIntensity = @(z_val, I) -alpha_proteins * I - TPA_prot * I.^2;
I_z = zeros(1, n_depth); % Initialize spatial intensity profile (time-independent)

% Solve for spatial intensity profile using ode45
[~, I_sol] = ode45(odeIntensity, z, I0);
I_z = I_sol'; % Store as row vector

% Create 2D intensity matrix (space x time)
I_z_t = zeros(n_depth, n_time);
for ti = 1:n_time
    % Scale depth profile by temporal profile
    I_z_t(:, ti) = I_z * (I_t(ti) / I0);
end

%% --- Calculate Photons per Chromophore per Pulse ---
tau_pulse = FWHM; % pulse duration [s]
Photon_per_chromo = zeros(1, n_depth);

% Calculate photon flux and photons per chromophore
gamma_flux = I_z * tau_pulse / E_photon; % [photons/m²]
Photon_per_chromo = gamma_flux ./ (n_prot * z); % [photons/chromophore]
% Fix first point to avoid division by zero
Photon_per_chromo(1) = Photon_per_chromo(2);

%% --- Population Dynamics Calculations ---
% Population matrices
n0_matrix = zeros(n_depth, n_time); % Ground state
n1_matrix = zeros(n_depth, n_time); % First excited state
n4_matrix = zeros(n_depth, n_time); % Higher excited state
n01_matrix = zeros(n_depth, n_time); % 1PA relaxed state
n02_matrix = zeros(n_depth, n_time); % 2PA relaxed state
n1_dark_matrix = zeros(n_depth, n_time); % Dark state 1
n2_dark_matrix = zeros(n_depth, n_time); % Dark state 2

% Initialize the ground state population at t=0 for all depths
n0_matrix(:,1) = n_prot;

% ODE solver options
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);

% Pulse end time index (useful for after-pulse calculations)
pulse_end_idx = find(t > FWHM*3, 1);
if isempty(pulse_end_idx)
    pulse_end_idx = n_time;
end

% Solve rate equations for each depth
for zi = 1:n_depth
    % Initial values at this depth
    n0 = n_prot;
    n1 = 0;
    n4 = 0;
    n01 = 0;
    n02 = 0;
    n1_dark = 0;
    n2_dark = 0;
    
    % Store initial values
    n0_matrix(zi,1) = n0;
    n1_matrix(zi,1) = n1;
    n4_matrix(zi,1) = n4;
    n01_matrix(zi,1) = n01;
    n02_matrix(zi,1) = n02;
    
    % Create time-dependent intensity function for this depth
    I_func = @(t_in) interp1(t, I_z_t(zi,:), t_in);
    
    % First solve during the pulse using the full rate equations
    odefun_during_pulse = @(t_in, n) rate_eqs_during_pulse(t_in, n, t, I_func, sigma_01, sigma_04, lifetime_1, lifetime_4);
    
    % Solve for the pulse duration
    [t_sol1, n_sol1] = ode45(odefun_during_pulse, t(1:pulse_end_idx), [n0; n1; n4; n01; n02;], options);
    
    % If we need to continue after the pulse
    if pulse_end_idx < n_time
        % Get final state from pulse calculation
        final_state = n_sol1(end,:)';
        
        % Define ODE function for after-pulse period
        odefun_after_pulse = @(t_in, n) rate_eqs_after_pulse(t_in, n, lifetime_1, lifetime_4);
        
        % Solve for after-pulse period
        [t_sol2, n_sol2] = ode45(odefun_after_pulse, t(pulse_end_idx:end), final_state, options);
        
        % Combine solutions
        t_sol = [t_sol1; t_sol2(2:end)];
        n_sol = [n_sol1; n_sol2(2:end,:)];
    else
        t_sol = t_sol1;
        n_sol = n_sol1;
    end
    
    % Interpolate results back to our original time grid
    for ni = 1:5
        temp = interp1(t_sol, n_sol(:,ni), t);
        
        % Store in appropriate matrix
        switch ni
            case 1
                n0_matrix(zi,:) = temp;
            case 2
                n1_matrix(zi,:) = temp;
            case 3
                n4_matrix(zi,:) = temp;
            case 4
                n01_matrix(zi,:) = temp;
            case 5
                n02_matrix(zi,:) = temp;
        end
    end
end

%% --- Calculate 2PA/1PA Ratio ---
epsilon = 1e-30; % Small value to avoid division by zero
final_time_index = size(n01_matrix, 2); % Get the last time index

% Extract the final populations at each depth
final_1PA = n01_matrix(:, final_time_index);
final_2PA = n02_matrix(:, final_time_index);

% Calculate the ratio of 2PA to 1PA for each depth at final time
ratio_2PA_to_1PA = final_2PA ./ (final_1PA + epsilon);

% Create the ratio matrix for all depths and times
ratio_matrix = n02_matrix ./ (n01_matrix + epsilon);

% Weighted average of the 2PA/1PA Ratio
weighted_average_2PA_to_1PA = trapz(z, ratio_2PA_to_1PA) ./ (z(end) - z(1));

%Time dependent weighted averages: 
weighted_avg_ratio_vs_time = zeros(1, length(t));
weighted_avg_1PA_vs_time = weighted_avg_ratio_vs_time;
weighted_avg_2PA_vs_time=weighted_avg_ratio_vs_time;
for ti = 1:length(t)
    ratio_at_t = ratio_matrix(:, ti);
    SPA_at_t =  n01_matrix(:, ti)./n_prot;
    TPA_at_t =  n02_matrix(:, ti)./n_prot;

    weighted_avg_ratio_vs_time(ti) = trapz(z, ratio_at_t) / (z(end) - z(1));
    weighted_avg_1PA_vs_time(ti)=trapz(z, SPA_at_t) ./ (z(end) - z(1));
    weighted_avg_2PA_vs_time(ti)=trapz(z, TPA_at_t) ./ (z(end) - z(1));
end

%% --- Visualizations ---
% Composite plot 
figure('Name', 'Composite: n01, n02, 2PA vs 1PA vs Depth', 'Color', 'white', 'Position', [100, 100, 800, 900]);
hold on
% Convert depth to micrometers
depth_um = z * 1e6;

% Subplot 1: n01 vs depth
subplot(3,1,1);
plot(t*1e15, weighted_avg_1PA_vs_time, 'b-', 'LineWidth', 2);
xlabel('Time (fs)');
ylabel('n_{01} / n_{prot}');
title(sprintf('I3^-: %.0f GW/cm^2, %.0f fs, %.2f J/cm^2', Peak_power_FWHM./10^9,FWHM.*10^15,Peak_power_FWHM.*FWHM),'Electron Population in 1PA pathway');
grid on;

% Subplot 2: n02 vs depth
subplot(3,1,2);
plot(t*1e15,weighted_avg_2PA_vs_time, 'r-', 'LineWidth', 2);
xlabel('Time (fs)');
ylabel('n_{02} / n_{prot}');
title('','Electron Population in 2PA pathway');
grid on;

% Subplot 3: 3D 2PA/1PA vs depth vs time
subplot(3,1,3);
imagesc(t*1e15, z*1e6, ratio_matrix);
xlabel('Time (fs)'); ylabel('Depth (μm)');
title('','2PA/1PA Ratio');
colorbar; colormap(jet);grid on;

%% Print statistics
fprintf('Results:\n');
fprintf('Photons per chromophore at surface: %.2f\n', Photon_per_chromo(1));
fprintf('Photons per chromophore at 5 μm: %.2f\n', interp1(z*1e6, Photon_per_chromo, 5, 'linear', 'extrap'));
fprintf('Photons per chromophore at end depth: %.2f\n', Photon_per_chromo(end));
fprintf('Weighted average of 2PA/1PA: %.4f\n', weighted_average_2PA_to_1PA);
fprintf('Time constants used:\n');
fprintf('  - Lifetime S1: %.2f fs\n', lifetime_1*1e15);
fprintf('  - Lifetime S2: %.2f fs\n', lifetime_4*1e15);

%% --- Rate Equation Functions ---
% During pulse rate equations
function dn_dt = rate_eqs_during_pulse(t_solve, n, t_grid, I_func, sigma_01, sigma_04, lifetime_1, lifetime_4)
    % Get current intensity
    I = I_func(t_solve);
    
    % Extract state populations
    n0 = n(1);
    n1 = n(2);
    n4 = n(3);
    n01 = n(4);
    n02 = n(5);
    
    % Compute rate equations based on the provided model
    N_e = 0; % No external excitation
    
    dn0_dt = -N_e - I * (n0 - n1) * sigma_01 - I^2 / 2 * (n0 - n4) * sigma_04;
    dn1_dt = I * (n0 - n1) * sigma_01 - I * (n1 - n4) * sigma_01- n1/lifetime_1;
    dn4_dt = I^2 / 2 * (n0 - n4) * sigma_04 +I * (n1 - n4) * sigma_01- n4/(lifetime_4);
    dn01_dt = n1/lifetime_1;
    dn02_dt = n4/lifetime_4;
    
    % Return as column vector
    dn_dt = [dn0_dt; dn1_dt; dn4_dt; dn01_dt; dn02_dt];
end

% After pulse rate equations
function dn_dt = rate_eqs_after_pulse(t_solve, n, lifetime_1, lifetime_4)
    % Extract state populations
    n0 = n(1);
    n1 = n(2);
    n4 = n(3);
    n01 = n(4);
    n02 = n(5);
    
    % Compute rate equations based on the provided model for after-pulse period
    dn0_dt = 0;
    dn1_dt = -n1/lifetime_1;
    dn4_dt = -n4/lifetime_4;
    dn01_dt = n1/lifetime_1;
    dn02_dt = n4/lifetime_4;
    
    % Return as column vector
    dn_dt = [dn0_dt; dn1_dt; dn4_dt; dn01_dt; dn02_dt;];

end
