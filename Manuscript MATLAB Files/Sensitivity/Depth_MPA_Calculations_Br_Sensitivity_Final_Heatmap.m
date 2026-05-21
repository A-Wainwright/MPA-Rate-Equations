%% ============================================================
%  Br MRE HEATMAP ANALYSIS
%  Power vs τ_Detune and Power vs σ14
%  Power sweep: logspace(1e9 → 1e12 W/cm^2)
%% ============================================================

clear; clc; close all;

%% --- Constants ---
c = 299792458;
h_bar = 1.054571817e-34;

wavelength = 530e-9;
w = 2*pi*c/wavelength;
E_photon = h_bar*w;

%% --- Material (Br) ---
M_prot_solution = 0.3333e-3;
NA = 6.022e23;
n_prot = M_prot_solution * NA * 1e3;

molar_extinction_Protein = 6e4;
ext_coeff = molar_extinction_Protein * M_prot_solution;
alpha_proteins = ext_coeff * 100;

Cross_section_1 = alpha_proteins / n_prot;

sigma01_base = Cross_section_1 / (h_bar*w);
sigma04_base = (290e-58) / (h_bar*w)^2;
sigma14_base = sigma01_base;

TPA_base = (290e-58)/(h_bar*w) * n_prot;

%% --- Lifetimes ---
tau1 = 500e-15;
tau4 = 200e-15;

%% --- Parameter Sweeps ---
tau_values = linspace(200e-15, 500e-15, 20);   % τ_detune
sigma14_factors = linspace(0.5, 1.5, 20);      % σ14 scaling

% 🔥 LOG POWER SWEEP (requested)
power_cm = logspace(9, 12, 20);  % W/cm^2
power = power_cm * 1e4;          % convert → W/m^2

%% --- Grid ---
FWHM = 200e-15;

t = linspace(-2*FWHM, 5e-12, 100);
z = linspace(0.1e-6, 10e-6, 200);

options = odeset('RelTol',1e-4,'AbsTol',1e-6);

%% ============================================================
% HEATMAP 1: Power vs τ_detune
%% ============================================================

heat_tau = zeros(length(power), length(tau_values));

for ip = 1:length(power)
    I0 = power(ip);

    % Temporal profile
    I_t = I0 * 2.^(-(2*(t/FWHM)).^2);

    for itau = 1:length(tau_values)

        tau_d = tau_values(itau);

        heat_tau(ip, itau) = run_case_heat( ...
            I0, I_t, ...
            1, 1, ... % σ fixed
            tau_d, ...
            sigma01_base, sigma04_base, sigma14_base, ...
            alpha_proteins, n_prot, tau1, tau4, ...
            TPA_base, t, z, options);

    end
end

%% ============================================================
% HEATMAP 2: Power vs σ14
%% ============================================================

heat_sigma14 = zeros(length(power), length(sigma14_factors));

for ip = 1:length(power)
    I0 = power(ip);

    I_t = I0 * 2.^(-(2*(t/FWHM)).^2);

    for is = 1:length(sigma14_factors)

        f14 = sigma14_factors(is);

        heat_sigma14(ip, is) = run_case_heat( ...
            I0, I_t, ...
            1, f14, ... % vary σ14
            tau1, ...
            sigma01_base, sigma04_base, sigma14_base, ...
            alpha_proteins, n_prot, tau1, tau4, ...
            TPA_base, t, z, options);

    end
end

%% ============================================================
% PLOTS
%% ============================================================

figure('Color','white','Name','Br Heatmaps');

subplot(1,2,1)
imagesc(tau_values*1e15, log10(power_cm), heat_tau);
set(gca,'YDir','normal');
xlabel('\tau_{Detune} (fs)', 'FontSize', 14);
ylabel('log_{10}(Power) [W/cm^2]', 'FontSize', 14);
title('2PA (%) vs Power & \tau_{Detune}');
cb=colorbar;
colormap(jet);
cb.Label.String = '2PA (%)';   % or whatever quantity you're plotting
cb.Label.FontSize = 14;

subplot(1,2,2)
imagesc(sigma14_factors, log10(power_cm), heat_sigma14);
set(gca,'YDir','normal');
xlabel('\sigma_{14} scaling', 'FontSize', 14);
ylabel('log_{10}(Power) [W/cm^2]', 'FontSize', 14);
title('2PA (%) vs Power & \sigma_{14}');
cb=colorbar;
colormap(jet);
cb.Label.String = '2PA (%)';   % or whatever quantity you're plotting
cb.Label.FontSize = 14;

%sgtitle('Bacteriorhodopsin Sensitivity Maps');

%% ============================================================
% CORE HEATMAP FUNCTION
%% ============================================================

function avg_percent = run_case_heat( ...
    I0, I_t, ...
    f04, f14, tau_dephase, ...
    sigma01_base, sigma04_base, sigma14_base, ...
    alpha, n_prot, tau1, tau4, ...
    TPA_base, t, z, options)

    sigma01 = sigma01_base;
    sigma04 = f04 * sigma04_base;
    sigma14 = f14 * sigma14_base;

    TPA = f04 * TPA_base;

    % Depth-dependent intensity
    odeI = @(z,I) -alpha*I - TPA*I.^2;
    [~,I_sol] = ode45(odeI, z, I0);
    I_z = I_sol';

    % Build I(z,t)
    I_z_t = zeros(numel(z), numel(t));
    for ti = 1:numel(t)
        I_z_t(:,ti) = I_z * (I_t(ti)/I0);
    end

    n01 = zeros(numel(z),1);
    n02 = zeros(numel(z),1);

    for zi = 1:numel(z)

        I_func = @(tq) interp1(t, I_z_t(zi,:), tq,'linear','extrap');

        odefun = @(tq,n) rate_eqs_full(tq,n,I_func,...
            sigma01,sigma04,sigma14,tau1,tau4,tau_dephase);

        [~,n_sol] = ode15s(odefun, t, [n_prot;0;0;0;0;0;0], options);

        n01(zi) = n_sol(end,4);
        n02(zi) = n_sol(end,5);
    end

    eps = 1e-30;
    percent = 100 * n02 ./ (n01 + n02 + eps);

    avg_percent = trapz(z, percent)/(z(end)-z(1));
end

%% ============================================================
% RATE EQUATIONS
%% ============================================================

function dn = rate_eqs_full(t,n,I_func,sigma01,sigma04,sigma14,tau1,tau4,tau_dephase)

    I = I_func(t);

    n0=n(1); n1=n(2); n4=n(3);
    n01=n(4); n02=n(5);
    n1d=n(6); n2d=n(7);

    dn0 = -I*(n0-n1)*sigma01 - (I^2/2)*(n0-n4)*sigma04;
    dn1 = I*(n0-n1)*sigma01 - I*(n1-n4)*sigma14 - n1/tau_dephase;
    dn4 = (I^2/2)*(n0-n4)*sigma04 + I*(n1-n4)*sigma14 - n4/(tau_dephase*0.7);

    dn01 = n1d/tau1;
    dn02 = n2d/tau4;

    dn1d = n1/tau_dephase - n1d/tau1;
    dn2d = n4/(tau_dephase*0.7) - n2d/tau4;

    dn = [dn0;dn1;dn4;dn01;dn02;dn1d;dn2d];
end