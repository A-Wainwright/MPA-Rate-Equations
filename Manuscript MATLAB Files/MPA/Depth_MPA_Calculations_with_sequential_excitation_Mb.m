%% Multi Rate Equation (MRE) avec dépendance spatiale I(z,t)
% Version corrigée avec intensité en profondeur appropriée
clear; clc; close all;

%% --- Constants and Parameters ---
c = 299792458; % speed of light in m/s
e = 1.60217663e-19; % charge of the electron in C
h_bar = 1.054571817e-34; % reduced Planck's constant in J/Hz

% Wavelength and angular frequency
wavelength = 530e-9; % [m]
w = 2 * pi * c / wavelength;
E_photon = h_bar * w; % photon energy [J]

% Protein density and cross sections
M_prot_solution = 6e-3; % [M]
NA = 6.022e23; % Avogadro's number
n_prot = M_prot_solution * NA * 1e3; % [1/m^3]
molar_extinction_Protein = 11600; % M-1 cm-1 at 530 nm

% Cross sections calculation
extinction_coefficient_proteins = molar_extinction_Protein * M_prot_solution; % [1/cm]
Cross_section_1 = extinction_coefficient_proteins * 100/n_prot; % [m^2]
alpha_proteins = extinction_coefficient_proteins * 100; % [1/m]

% Experimental cross sections
Cross_section_0_1_per_photon = Cross_section_1; % [m^2]
Cross_section_0_4_per_photon = 100e-58; % [m^4*s/photon]
Cross_section_1_4_per_photon = Cross_section_0_1_per_photon; % [m^2]

% Convert to intensity-based cross sections
sigma_01 = Cross_section_0_1_per_photon / (h_bar * w); % [m^2 / (W/m^2)]
sigma_04 = Cross_section_0_4_per_photon / (h_bar * w)^2; % [m^4 / (W/m^2)^2]
sigma_14 = Cross_section_1_4_per_photon / (h_bar * w); % [m^2 / (W/m^2)]

% Two-photon absorptivity [1/(m·W)]
TPA_prot = Cross_section_0_4_per_photon / (h_bar * w) * n_prot;

% Lifetimes
tau_1 = 70e-15; % [s]
tau_2 = 50e-15; % [s]
dephasing = tau_1;

% Laser pulse parameters
FWHM = 80e-15; % [s]
Peak_power_FWHM = 250e9; % [W/cm^2]
I0_cm = Peak_power_FWHM; % [W/cm^2]
I0 = I0_cm * 1e4; % [W/m^2]

% Time and depth grid
t_max = 2e-12; % Maximum time [s]
t_pre = -2 * FWHM; % Start time [s]
n_time = 500; % Number of time points
z_max = 10e-6; % Maximum depth [m]
z_min = 0.1e-6; % Minimum depth [m]
n_depth = 300; % Number of depth points

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

% Initialize the ground state population at t=0 for all depths
n0_matrix(:,1) = n_prot;

% ODE solver options
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);

% Solve rate equations for each depth
for zi = 1:n_depth
    % Initial values at this depth
    n0 = n_prot;
    n1 = 0;
    n4 = 0;
    n01 = 0;
    n02 = 0;
    
    % Store initial values
    n0_matrix(zi,1) = n0;
    n1_matrix(zi,1) = n1;
    n4_matrix(zi,1) = n4;
    n01_matrix(zi,1) = n01;
    n02_matrix(zi,1) = n02;
    
    % Create time-dependent intensity function for this depth
    I_func = @(t_in) interp1(t, I_z_t(zi,:), t_in);
    
    % Define ODE function with time-dependent intensity
    odefun = @(t_in, n) rate_eqs(n, I_func(t_in), sigma_01, sigma_14, sigma_04, tau_1, tau_2, dephasing);
    
    % Solve for all times at once
    [t_sol, n_sol] = ode15s(odefun, t, [n0; n1; n4; n01; n02], options);
    
    % Store results
    n0_matrix(zi,:) = n_sol(:,1)';
    n1_matrix(zi,:) = n_sol(:,2)';
    n4_matrix(zi,:) = n_sol(:,3)';
    n01_matrix(zi,:) = n_sol(:,4)';
    n02_matrix(zi,:) = n_sol(:,5)';
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
weighted_average_2PA_to_1PA=trapz(z, ratio_2PA_to_1PA) ./ (z(end) - z(1));

%% --- Visualizations ---
% Plot intensity profile vs depth
figure('Name', 'Intensity vs Depth', 'Color', 'white');
plot(z*1e6, I_z*1e-4*1e-9, 'LineWidth', 1.5);
xlabel('Depth (μm)');
ylabel('Intensity (GW/cm^2)');
title('Intensity Profile vs Depth');
grid on;

% Plot photons per chromophore vs depth
figure('Name', 'Photons per Chromophore vs Depth', 'Color', 'white');
semilogy(z*1e6, Photon_per_chromo, 'LineWidth', 1.5);
xlabel('Depth (μm)');
ylabel('Photons per chromophore per pulse');
title('Photon Flux per Chromophore vs. Depth');
grid on;

% Time-resolved depth profiles
figure('Name', 'Space-Time Population Profiles', 'Color', 'white');
subplot(2,2,1)
imagesc(t*1e15, z*1e6, n01_matrix/n_prot);
xlabel('Time (fs)'); ylabel('Depth (μm)');
title('1PA Normalized Population'); 
colorbar; colormap(jet);

subplot(2,2,2)
imagesc(t*1e15, z*1e6, n02_matrix/n_prot);
xlabel('Time (fs)'); ylabel('Depth (μm)');
title('2PA Normalized Population'); 
colorbar; colormap(jet);

subplot(2,2,3)
imagesc(t*1e15, z*1e6, n01_matrix./max(n02_matrix(:)));
xlabel('Time (fs)'); ylabel('Depth (μm)');
title('1PA Population (Normalized to max 2PA)'); 
colorbar; colormap(jet);

subplot(2,2,4)
imagesc(t*1e15, z*1e6, I_z_t/I0);
xlabel('Time (fs)'); ylabel('Depth (μm)');
title('Intensity I(z,t)/I_0'); 
colorbar; colormap(hot);

% 2PA/1PA ratio vs depth
figure('Name', '2PA/1PA Ratio vs Depth', 'Color', 'white');
semilogy(z*1e6, ratio_2PA_to_1PA, 'b-o', 'LineWidth', 1.5);
yline(weighted_average_2PA_to_1PA,'-','Average');

xlabel('Depth (μm)');
ylabel('2PA/1PA Ratio');
title(sprintf('MbCO: %.2f GW/cm^2', Peak_power_FWHM./10^9));
grid on;

% combined figure
figure('Name', 'Space-Time Population Profiles', 'Color', 'white');
subplot(1,2,1)
imagesc(t*1e15, z*1e6, ratio_matrix);
xlabel('Time (fs)'); ylabel('Depth (μm)');
title(sprintf('MbCO: 2PA/1PA Ratio for %.2f GW/cm^2', Peak_power_FWHM./10^9));
colorbar; colormap(jet);
subplot(1,2,2)
semilogy(z*1e6, ratio_2PA_to_1PA, 'b-o', 'LineWidth', 1.5);
yline(weighted_average_2PA_to_1PA,'-','Average');
xlabel('Depth (μm)');
ylabel('2PA/1PA Ratio');
title(sprintf('MbCO:10 ps 2PA/1PA Ratio for  %.2f GW/cm^2', Peak_power_FWHM./10^9));
grid on;


% Print statistics
fprintf('Results:\n');
fprintf('Photons per chromophore at surface: %.2f\n', Photon_per_chromo(1));
fprintf('Photons per chromophore at 5 μm: %.2f\n', interp1(z*1e6, Photon_per_chromo, 5));
fprintf('Photons per chromophore at 10 μm: %.2f\n', Photon_per_chromo(end));
fprintf('weighted average of 2PA/1PA %.4f\n', weighted_average_2PA_to_1PA);

%% --- Rate Equation Function ---
function dn_dt = rate_eqs(n, I, sigma_01, sigma_14, sigma_04, tau_1, tau_4, dephasing)
    n0 = n(1);
    n1 = n(2);
    n4 = n(3);
    n01 = n(4);
    n02 = n(5);

    % Pre-calculate differences to avoid redundant calculations
    n0_n1_diff = n0 - n1;
    n0_n4_diff = n0 - n4;
    n1_n4_diff = n1 - n4;
    
    % Rate equations
    dn0 = -I * n0_n1_diff * sigma_01 - I^2/2 * n0_n4_diff * sigma_04;
    dn1 = I * n0_n1_diff * sigma_01 - I * n1_n4_diff * sigma_14 - n1/dephasing;
    dn4 = I^2/2 * n0_n4_diff * sigma_04 + I * n1_n4_diff * sigma_14 - n4/tau_4;
    dn01 = n1 / tau_1;
    dn02 = n4 / tau_4;

    dn_dt = [dn0; dn1; dn4; dn01; dn02];
end