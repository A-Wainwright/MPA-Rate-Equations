%% Rate Equation modeling for MB
% 4-panel sensitivity study + composite figure

clear; clc; close all;

%% --- Constants and Parameters ---
c     = 299792458;
h_bar = 1.054571817e-34;

wavelength = 530e-9;
w = 2*pi*c/wavelength;
E_photon = h_bar*w;

M_prot_solution = 6e-3;
NA = 6.022e23;
n_prot = M_prot_solution * NA * 1e3;

molar_extinction_Protein = 11600;

extinction_coefficient_proteins = molar_extinction_Protein * M_prot_solution;
alpha_proteins  = extinction_coefficient_proteins * 100;

Cross_section_1 = extinction_coefficient_proteins * 100 / n_prot;

Cross_section_0_1_per_photon_base = Cross_section_1;
Cross_section_0_4_per_photon_base = 100e-58;
Cross_section_1_4_per_photon_base = Cross_section_0_1_per_photon_base;

%% --- Time constants ---
tau_1_base = 70e-15;
tau_4_base = 50e-15;
dephasing  = tau_1_base;

%% --- Laser ---
FWHM = 80e-15;
Peak_power_FWHM = 100e9;
I0 = Peak_power_FWHM * 1e4;

%% --- Grid ---
t_max  = 2e-12;
t_pre  = -2 * FWHM;
n_time = 500;

z_max = 10e-6;
z_min = 0.1e-6;
n_depth = 300;

t = linspace(t_pre, t_max, n_time);
z = linspace(z_min, z_max, n_depth);

I_t = I0 * 2.^(-(2*(t/FWHM)).^2);

options = odeset('RelTol',1e-4,'AbsTol',1e-6);

%% ============================================================
% SENSITIVITY SETTINGS
%% ============================================================

factors = [0.5 0.9 1.0 1.1 1.5];
nf = numel(factors);

tau_factors = [0.5 0.75 1.0 1.25 1.5];
ntau = numel(tau_factors);

%% --- STORAGE ---
percent_diag   = zeros(n_depth,nf);
percent_f0neq1 = zeros(n_depth,nf);
percent_f1neq1 = zeros(n_depth,nf);
percent_tau    = zeros(n_depth,ntau);

%% ============================================================
% RUN: σ0n & σ1n TOGETHER
%% ============================================================
for k = 1:nf

    f0n = factors(k);
    f1n = factors(k);

    [p,~] = run_case( ...
        f0n,f1n,1,1, ...
        Cross_section_0_1_per_photon_base, ...
        Cross_section_0_4_per_photon_base, ...
        Cross_section_1_4_per_photon_base, ...
        alpha_proteins,n_prot,h_bar,w,E_photon,...
        tau_1_base,tau_4_base,I0,I_t,t,z,FWHM,options);

    percent_diag(:,k) = p;
end

%% ============================================================
% RUN: σ0n sweep
%% ============================================================
for k = 1:nf

    f0n = factors(k);
    f1n = 1.0;

    [p,~] = run_case( ...
        f0n,f1n,1,0, ...
        Cross_section_0_1_per_photon_base, ...
        Cross_section_0_4_per_photon_base, ...
        Cross_section_1_4_per_photon_base, ...
        alpha_proteins,n_prot,h_bar,w,E_photon,...
        tau_1_base,tau_4_base,I0,I_t,t,z,FWHM,options);

    percent_f0neq1(:,k) = p;
end

%% ============================================================
% RUN: σ1n sweep
%% ============================================================
for k = 1:nf

    f0n = 1.0;
    f1n = factors(k);

    [p,~] = run_case( ...
        f0n,f1n,0,1, ...
        Cross_section_0_1_per_photon_base, ...
        Cross_section_0_4_per_photon_base, ...
        Cross_section_1_4_per_photon_base, ...
        alpha_proteins,n_prot,h_bar,w,E_photon,...
        tau_1_base,tau_4_base,I0,I_t,t,z,FWHM,options);

    percent_f1neq1(:,k) = p;
end

%% ============================================================
% COMPOSITE FIGURE
%% ============================================================

figure('Color','white');

cols = lines(nf);

all_data = [percent_diag(:); percent_f0neq1(:); percent_f1neq1(:); percent_tau(:)];

ymin = max(min(all_data(all_data>0))*0.9,1e-12);
ymax = max(all_data)*1.1;

%% --- (1,1) σ0n & σ1n together ---
subplot(3,1,1); hold on;
for k = 1:nf
    semilogy(z*1e6,max(percent_diag(:,k),1e-12), ...
        'Color',cols(k,:), ...
        'LineWidth',1.6, ...
        'DisplayName',sprintf('%.1fx',factors(k)));
end
title('\sigma_{0n} & \sigma_{1n}');
xlabel('Depth (μm)'); ylabel('2PA (%)');
ylim([ymin ymax]); grid on; legend;

%% --- (1,2) σ0n ---
subplot(3,1,2); hold on;
for k = 1:nf
    semilogy(z*1e6,max(percent_f0neq1(:,k),1e-12), ...
        'Color',cols(k,:), ...
        'LineWidth',1.6, ...
        'DisplayName',sprintf('σ_{0n}=%.1fx',factors(k)));
end
title('\sigma_{0n} sensitivity');
xlabel('Depth (μm)'); ylabel('2PA (%)');
ylim([ymin ymax]); grid on; legend;

%% --- (2,1) σ1n ---
subplot(3,1,3); hold on;
for k = 1:nf
    semilogy(z*1e6,max(percent_f1neq1(:,k),1e-12), ...
        'Color',cols(k,:), ...
        'LineWidth',1.6, ...
        'DisplayName',sprintf('σ_{1n}=%.1fx',factors(k)));
end
title('\sigma_{1n} sensitivity');
xlabel('Depth (μm)'); ylabel('2PA (%)');
ylim([ymin ymax]); grid on; legend;

%% =======================================
%%% --- Global title ---
sgtitle(sprintf('MRE Sensitivity Analysis (Peak Power = %.2f GW/cm^2)', ...
    I0/(1e9*10^4)));
%% =====================

%% =====================
% FUNCTION
%% ============================================================

function [percent_2PA_depth, weighted_avg] = run_case( ...
    f0n,f1n,~,~, ...
    C01,C0n,C1n, ...
    alpha,n_prot,h_bar,w,E, ...
    tau1,tau4,I0,I_t,t,z,FWHM,options)

sigma_01 = C01/(h_bar*w);
sigma_0n = f0n*C0n/(h_bar*w)^2;
sigma_1n = f1n*C1n/(h_bar*w);

TPA = C0n/(h_bar*w)*n_prot;

odeI = @(z,I) -alpha*I - TPA*I.^2;
[~,Isol] = ode45(odeI,z,I0);
I_z = Isol.';

I_zt = zeros(numel(z),numel(t));
scale = I_t/I0;

for i=1:numel(t)
    I_zt(:,i)=I_z.'*scale(i);
end

n01=zeros(numel(z),numel(t));
n02=zeros(numel(z),numel(t));

for zi=1:numel(z)

    If = @(tt) interp1(t,I_zt(zi,:),tt);

    ode = @(tt,n) rate_eq(n,If(tt),sigma_01,sigma_1n,sigma_0n,tau1,tau4);

    [~,sol]=ode15s(ode,t,[n_prot;0;0;0;0],options);

    n01(zi,:)=sol(:,4).';
    n02(zi,:)=sol(:,5).';
end

eps0=1e-30;
p = 100*n02(:,end)./(n01(:,end)+n02(:,end)+eps0);

weighted_avg = trapz(z,p)/(z(end)-z(1));

percent_2PA_depth = p;
end

function dn = rate_eq(n,I,s01,s1n,s0n,t1,t4)

n0=n(1); n1=n(2); n4=n(3);

dn0 = -I*(n0-n1)*s01 - 0.5*I^2*(n0-n4)*s0n;
dn1 =  I*(n0-n1)*s01 - I*(n1-n4)*s1n - n1/t1;
dn4 = 0.5*I^2*(n0-n4)*s0n + I*(n1-n4)*s1n - n4/t4;

dn01 = n1/t1;
dn02 = n4/t4;

dn = [dn0;dn1;dn4;dn01;dn02];
end