%% Protein SRE Model
%By Alexander Wainwright 
clear;clc;close all;

%% Notes 
%Need to fix the attenuation caused by the protein vs the liberated
%electrons! line 129 and after

%% Constants
c = 299792458; %speed of light in m./s
e= 1.60217663.*10.^-19; %charge of the electron in C
m = 9.11.*10.^-31; %mass of electron
m_c = m./2;%mass of electron in conduction band
%M_h20= 4.692.*10.^-23;%[Kg]
mu0=4.*pi.*10.^-7;
E_0 = 8.854.*10.^-12; %permitivity of free sapce in F./m.
E_0_ev =55.26349406.*10.^6; %permittivity of free space in e.^2⋅eV.^−1⋅m−.^1;
h_bar = 1.054571817.*10.^-34; %Plank's reduced Constant in J./Hz
h_bar_ev = 6.5821.*10.^-16; %Plank's reduced Constant in Ev.*s
n_prot_cm=1.81.*10.^15;
n_prot=n_prot_cm.*10.^6;%1./m.^3, 1 electron per protein

molar_extinction_Protein=63000; % M-1 cm-1 at 560 nm

n_bound_H20 = 10.^22.*10.^6;%1./m.^3, Stauration of bound water
n_bound_prot=(4044-1).*n_prot;%1./m.^3, 1 bond per atom per protein minus 1.

n_ini_prot = n_prot;%1./m.^3, 1 electron per protein
rho_ini_max_H20 = 10.^19.*10.^6;%m.^-3 intermediate level stauartion of water
%% Variables 
wavelength_1 =530.*10.^-9; %m
m_prime = 1./2.*m; %effective mass 
E_int_ev_protein=2.21; %intermediate gap in eV from absorption spectrum at 560 nm
E_gap_eV_protein = 7; %band gap in eV for S5 excitation
E_gap_eV_H20=9.5; %band gap in eV
E_int_H20_eV=6.4; %eV

E_gap_protein = E_gap_eV_protein.*1.60218e-19; %band gap in J
E_int_protein= E_int_ev_protein.*1.60218e-19; %intermeidate band gap in J
E_gap_H20=E_gap_eV_H20.*1.60218e-19; %band gap in J
E_int_H20=E_int_H20_eV.*1.60218e-19; %band gap in J

% Laser intensities
Peak_power_FWHM=4225*10.^9;
I_max_cm_1 = 2*sqrt(log(2))/sqrt(pi)*Peak_power_FWHM; %incident laser intensity W./cm.^2
I_max_1=I_max_cm_1.*10.^4; %incident laser intensity W./m.^2

%Laser pulse duration
T_laser_1 = 100.*10.^-15; %s 

%PPC
Cross_section_protein = 8.14*10^-21;
E_p=h_bar*2.*pi.*c./wavelength_1;
F= Peak_power_FWHM*(10^4)*T_laser_1;
PPC=F/E_p*Cross_section_protein;


%Material constants
t_coll = 1.*10.^-15; %s
% laser-induced conduction band electrons in water. By Vogel
% Source for ns ablation https:././www.ncbi.nlm.nih.gov./pmc./articles./PMC6422691./#

%% Space and time dependent intensities
length_measure=2.*T_laser_1;
%t=linspace(-lenght_measure,lenght_measure,500);

L = length_measure;  % rename for brevity
t_first  = linspace(-L, 0, 300);       % includes 0 as the last point
t_second = linspace(0,  L, 101);       % includes 0 as the first point
t = [t_first, t_second(2:end)];        % drop the duplicate 0

% Time dependent intensity with gausien shape
t_0_1=0;
t_0_2=T_laser_1;
I_1= I_max_1.*2.^(-1.*(2.*(t-t_0_1./2)./T_laser_1).^2); %time dependent intensity profile.

%% index of refraction calculation for water
% Best guess
n_0_protein= 1.5;

% index of refraction calculation for water
% from:M. N. Polyanskiy. Refractiveindex.info database of optical constants. Sci. Data 11, 94 (2024)
%https:././doi.org./10.1038./s41597-023-02898-2
wavelength_um_1=wavelength_1.*10.^6; %um
n_0_H20=sqrt(1+5.666959820E-1./(1-5.084151894E-3./wavelength_um_1.^2)+1.731900098E-1./(1-1.818488474E-2./wavelength_um_1.^2)+2.095951857E-2./(1-2.625439472E-2./wavelength_um_1.^2)+1.125228406E-1./(1-1.073842352E1./wavelength_um_1.^2));

%% Rate Equation Canstants for wavelength 1
w = 2.*pi.*c./wavelength_1; %angular frequency of incident light 

F_Protein=sqrt(2.*I_1./(n_0_protein.*E_0.*c));% Electric Feild Intensity in protein in V./cm
gama_gap_protein = w./e.*sqrt(m_prime.*c.*n_0_protein.*E_0.*E_gap_protein./I_1); %Kedish peramiter
gama_int_protein = w./e.*sqrt(m_prime.*c.*n_0_protein.*E_0.*E_int_protein./I_1); %Kedish peramiter
%Cross_section_protein= t_coll./(w.^2.*t_coll.^2+1).*e.^2./(c.*n_0_protein.*E_0.*m_c); %Single photon colission cross section in m^2
W_1pt_protein = Cross_section_protein.*I_1./(h_bar.*w);

F_H20=sqrt(2.*I_1./(n_0_H20.*E_0.*c));% Electric Feild Intensity in H20 in V./cm
gama_gap_H20 = w./e.*sqrt(m_prime.*c.*n_0_H20.*E_0.*E_gap_H20./I_1); %Kedish peramiter
gama_int_H20 = w./e.*sqrt(m_prime.*c.*n_0_H20.*E_0.*E_int_H20./I_1); %Kedish peramiter
Cross_section_H20= t_coll./(w.^2.*t_coll.^2+1).*e.^2./(c.*n_0_H20.*E_0.*m_c); %Single photon colission cross section in m
W_1pt_H20 = Cross_section_H20.*I_1./(h_bar.*w);

%% Delta tilda (solving for the mixed case of gamma greater and less than one)
[diff_n_sfi_protein,t_rho_init_protein,delta_tilda_ev_gap_protein] = delta_tilda(gama_gap_protein,gama_int_protein,E_gap_protein,E_int_H20,m,w,e,F_Protein,t,n_bound_prot);
[diff_n_sfi_H20,t_rho_init_H20,delta_tilda_ev_gap_H20] = delta_tilda(gama_gap_H20,gama_int_H20,E_gap_H20,E_int_H20,m,w,e,F_H20,t,rho_ini_max_H20);

%Plotting the SFI components
% figure(1)
% semilogy(t_rho_init_protein*10^15,diff_n_sfi_protein,'k:',LineWidth=1);%1./cm.^3
% hold on
% semilogy(t_rho_init_H20*10^15,diff_n_sfi_H20,"k--",LineWidth=1);%1./cm.^3
% 
% title("Rate of SFI")
% ylabel('SFI (1/s*cm^{-3})')
% xlabel('Time (fs)')
% legend({'SFI of Proteins','SFI of Water'}, 'Location','southeast')

%% Rate equation solvers
%total number of electrons
delta_tilda_gap_protein=delta_tilda_ev_gap_protein.*1.60218e-19;
n_AI_Asym_Protein = log(2).*t_coll./(w.^2.*t_coll.^2+1).*((e.^2.*I_1./(c.*n_0_protein.*E_0.*m_c.*3./2.*(delta_tilda_gap_protein)))); %Estimation from later papers: Wavelength dependence of femtosecond ... Result matches nicely

delta_tilda_gap_H20=delta_tilda_ev_gap_H20.*1.60218e-19;
n_AI_Asym_H20 = log(2).*t_coll./(w.^2.*t_coll.^2+1).*((e.^2.*I_1./(c.*n_0_H20.*E_0.*m_c.*3./2.*(delta_tilda_gap_H20)))); %Estimation from later papers: Wavelength dependence of femtosecond ... Result matches nicely

%SRE(1) H20
[t_n_total,n_total] = ode45(@(t_solve,n) SRE1(t_solve,n,t,diff_n_sfi_protein,n_AI_Asym_Protein,n_bound_prot,diff_n_sfi_H20,n_AI_Asym_H20,n_bound_H20),t,[0,0]);
n_total_cm_3=(n_total(:,1)+n_total(:,2)).*10^-6;

%SFI electrons Proteins
[t_SFI_protein,n_SFI_protein]=ode45(@(t_sln,n) interp1(t_rho_init_protein,diff_n_sfi_protein,t_sln),t,0);
n_SFI_protein_cm_3=n_SFI_protein.*10^-6;

%SFI electrons Water
[t_SFI_H20,n_SFI_H20]=ode45(@(t_sln,n) interp1(t_rho_init_H20,diff_n_sfi_H20,t_sln),t,0);
n_SFI_H20_cm_3=n_SFI_H20.*10^-6;

%AI electrons
n_SFI_cm_3_interp=interp1(t_SFI_protein,n_SFI_protein_cm_3,t_n_total)+interp1(t_SFI_H20,n_SFI_H20_cm_3,t_n_total);
n_AI_cm_3=n_total_cm_3-n_SFI_cm_3_interp;

%% Intensity attenuation through the crystal of interest: this part is not working
% M_protein=n_prot_cm*10^3/(6.02*10^23);
% extinction_coefficient_electrons_p = Cross_section_protein.*n_total(:,1); %1/m
% extinction_coefficient_electrons_H20 = Cross_section_H20.*n_total(:,2);%1/m
% extinction_coefficient_proteins = molar_extinction_Protein.*M_protein.*100;%1/m
% z_1 = linspace(0,6*10^-6,25); %thickness of the sample in m
% z_2=linspace(100*10^-9,6*10^-6,75); %thickness of the sample in m
% z=[z_1,z_2];
% 
% surface_attenuation=extinction_coefficient_electrons_p+extinction_coefficient_electrons_H20+extinction_coefficient_proteins; %(1/m)
% crystal_attenuation = extinction_coefficient_proteins; %(1/M 1/cm)
% 
% [Time,depth]=meshgrid(t,z);
% 
% I_attenuation_surface=I_1.*exp(-(surface_attenuation).*(z_1));
% I_crystal=I_attenuation_surface(:,end).*exp(-(crystal_attenuation).*(z_2));
% I_z = [I_attenuation_surface,I_crystal];

%% Plots 
figure()
yyaxis left
ax = gca; % Get current axes
ax.FontSize = 16; % Set font size

semilogy(t_n_total*10^15,n_total_cm_3,'k:',LineWidth=3);%1./cm.^3
hold on
semilogy(t_SFI_protein*10^15,n_AI_cm_3,"k--",LineWidth=3);%1./cm.^3
semilogy(t_SFI_H20*10^15,n_SFI_H20_cm_3,LineWidth=3);%1./cm.^3
semilogy(t_SFI_protein*10^15,n_SFI_protein_cm_3,LineWidth=3);%1./cm.^3
ylim([10^12,10^22])
xlim([-100,100])
ylabel('Electron density(cm^{-3})','FontName', 'Arial', 'FontSize', 12)

yyaxis right
semilogy(t*10^15,I_1,LineWidth=3)
ylim([10^15,10^17])

ylabel('Laser Intensity (W/cm^{2})','FontName', 'Arial', 'FontSize', 12)
title(sprintf('Ionized Electron Density (FWHM = %.0f fs, %.0f PPC)', ...
    T_laser_1*10^15, PPC),'FontName', 'Arial', 'FontSize', 12);

xlabel('Time (fs)','FontName', 'Arial', 'FontSize', 12)
legend({'Total Ionization','AI','SFI of Water','SFI of Proteins'}, 'Location','southeast','FontName', 'Arial', 'FontSize', 12)

%% Compositve figure with intensity
figure
tlo = tiledlayout(4,1, 'TileSpacing','compact', 'Padding','compact');

% ================= TOP (3 rows): Electron density =================
axTop = nexttile([3 1]);  % span 3 rows, 1 column
semilogy(axTop, t_n_total*1e15, n_total_cm_3, 'k:', 'LineWidth', 3); hold(axTop, 'on')
semilogy(axTop, t_SFI_protein*1e15, n_AI_cm_3, 'k--', 'LineWidth', 3);
semilogy(axTop, t_SFI_H20*1e15, n_SFI_H20_cm_3, 'LineWidth', 3);
semilogy(axTop, t_SFI_protein*1e15, n_SFI_protein_cm_3, 'LineWidth', 3);
grid(axTop, 'on')
set(axTop, 'YLim', [1e12 1e22], 'XLim', [-100 100], 'FontSize', 12)
ylabel(axTop, 'Electron density (cm^{-3})', 'FontName','Arial', 'FontSize', 12)
%title(axTop, sprintf('Ionized Electron Density (FWHM = %.0f fs, %.0f PPC)', ...
%    T_laser_1*1e15, PPC), 'FontName','Arial', 'FontSize', 12);
legend(axTop, {'Total Ionization','AI','SFI of Water','SFI of Proteins'}, ...
       'Location','southeast', 'FontName','Arial', 'FontSize', 12)

% ================= BOTTOM (1 row): Laser intensity =================
axBottom = nexttile;      % occupies the remaining 1 row
plot(axBottom, t*1e15, I_1, 'LineWidth', 3);
grid(axBottom, 'on')
set(axBottom, 'YLim', [1e15 1e17], 'XLim', [-100 100], 'FontSize', 12)
ylabel(axBottom, 'Intensity (W/cm^{2})', 'FontName','Arial', 'FontSize', 12)

% Shared X label for the whole figure
xlabel(tlo, 'Time (fs)', 'FontName','Arial', 'FontSize', 12);

% Keep x-axes aligned on zoom/pan
linkaxes([axTop, axBottom], 'x');
% figure(3)
% semilogy(t_SFI_H20*10^15,n_SFI_H20_cm_3,LineWidth=1);%1./cm.^3
% hold on
% semilogy(t_SFI_protein*10^15,n_SFI_protein_cm_3,LineWidth=1);%1./cm.^3
% 
% title("Single rate equation")
% ylabel('Electron density(cm^{-3})')
% xlabel('Time (fs)')
% legend({'SFI of Water','SFI of Proteins'}, 'Location','southeast')
% 
% figure(4)
% semilogy(t_SFI_H20*10^15,diff_n_sfi_H20,LineWidth=1);%1./cm.^3
% hold on
% semilogy(t_SFI_protein*10^15,diff_n_sfi_protein,LineWidth=1);%1./cm.^3
% 
% title("Single rate equation")
% ylabel('SFI Rate ( 1/s cm^{-3})')
% xlabel('Time (fs)')
% legend({'SFI of Water','SFI of Proteins'}, 'Location','southeast')

% figure(5)
% semilogy(t_SFI_H20*10^15,diff_n_sfi_H20,LineWidth=1);%1./cm.^3
% title("Single rate equation")
% ylabel('SFI Rate ( 1/s cm^{-3})')
% xlabel('Time (fs)')
% legend({'SFI of Water','SFI of Proteins'}, 'Location','southeast')
% 
% figure(6)
% semilogy(t,extinction_coefficient_electrons_p+extinction_coefficient_electrons_H20)
% hold on
% semilogy(t,extinction_coefficient_proteins.*ones(size(t)))
% semilogy(t,surface_attenuation)
% ylabel('Attenuation (cm^{-1})')
% xlabel('Time (fs)')
% 
% figure(7)
% surf(Time, depth, I_attenuation_surface);
% xlabel('Time (s)');
% ylabel('Depth (m)');
% zlabel('Intensity');
% title('Intensity Attenuation Over Time and Depth');
% shading interp; 
% colorbar;

%% Helper functions
% Full keldish formula for photoexication
function dn_dt = SRE1(t_solve, n, t, diff_n_sfi_protein, n_AI_Asym_Protein, n_bound_prot, ...
                      diff_n_sfi_H20, n_AI_Asym_H20, n_bound_H20)
                  
    % Define recombination rate (set to 0 as it is considered negligible)
    recombination_rate = 0; 

    % Interpolate external inputs at the current time step
    diff_n_sfi_protein_interp = interp1(t, diff_n_sfi_protein, t_solve, 'linear', 'extrap');   
    n_AI_Asym_Protein_interp = interp1(t, n_AI_Asym_Protein, t_solve, 'linear', 'extrap'); 

    diff_n_sfi_H20_interp = interp1(t, diff_n_sfi_H20, t_solve, 'linear', 'extrap');   
    n_AI_Asym_H20_interp = interp1(t, n_AI_Asym_H20, t_solve, 'linear', 'extrap');

    % Ensure n has two components
    n1 = n(1);
    n2 = n(2);

    % Compute rate equations
    dn1_dt = diff_n_sfi_protein_interp + (n_AI_Asym_Protein_interp) * (n1 + n2) * ...
             ((n_bound_prot - n1) / (n_bound_prot + n_bound_H20)) - ...
             recombination_rate * (n1 + n2)^2;

    dn2_dt = diff_n_sfi_H20_interp + n_AI_Asym_H20_interp * (n1 + n2) * ...
             ((n_bound_H20 - n2) / (n_bound_prot + n_bound_H20)) - ...
             recombination_rate * (n1 + n2)^2;

    % Return as a column vector
    dn_dt = [dn1_dt; dn2_dt];

end

function [q]=Q(gama,x)
% Equation defined in Mechanisms of femtosecond laser nanosurgery of cells
% and tissues page 6.

    [K_gama,E_gama]=ellipke(gama./(sqrt(1+gama.^2)));
    [K,E] = ellipke(1./(sqrt(1+gama.^2)));
    %top value for the infinite sum estimation
    top = 500;
    x_val = linspace(0,top,top+1);
    %Definition of the function that will be summed
    funct_Q=@(l) exp(-pi.*l.*(K_gama-E_gama)./E).*dawson(sqrt(pi.*(2.*floor(x+1)-2.*x+l)./(2.*K.*E)));

    Q_sum_val =0;
    %Summation process
    for k = x_val
        Q_sum_val = Q_sum_val+funct_Q(k);
    end 

    %Value being returned
    q=sqrt(pi./(2.*K)).*Q_sum_val;
end

function [n_gamma_factor_constant]=n_sfi(delta_tilda,gama,w) %E_feild in eV
%E_gap in J (not eV)
%Returns the rate equation for the excitation of electrons acroos the band
%gap of intrest as a result of the strong feild ionization.
% Quations from: Mechanisms of femtosecond laser nanosurgery of cells and tissues page 6
    m= 9.11.*10.^-31; %[kg]
    h_bar = 1.054571817.*10.^-34;
    
    [K_1,E_1]=ellipke(gama./(sqrt(1+gama.^2)));
    E_2 = ellipticE(1./(sqrt(1+gama.^2))); %elliptic integral of [the first, the second] kind
    x=delta_tilda./(h_bar.*w); %Number of photons needed to exceed the band gap
    Q_val=Q(gama,x);

    %Full Keldysh model
    n_gamma_factor_constant = 2.*w./(9.*pi).*((sqrt(1+gama.^2)./gama).*m.*w./h_bar).^(3./2).*Q_val.*exp(-pi.*floor(x+1).*(K_1-E_1)./E_2);
end

function dndt = dn_sft_solver(t,n,t_sln,n_sfi_E_ini,rho_ini_max)
% Returns the rate equations which need to be solved. Rate equations 
%are defined in: Wavelength dependence of femtosecond laser-induced 
%breakdown in water and implications for laser surgery.
    n_sfi_E_ini= interp1(t_sln,n_sfi_E_ini,t);   
    dndt=n_sfi_E_ini.*(1-n/rho_ini_max);
end

function [diff_n_sfi,t_rho_init,delta_tilda_ev_gap] = delta_tilda(gama_gap,gama_int,E_gap,E_int,m,w,e,F,t,n_sat)
    %function solving for the change in the band gat (Delta tilda) for a given
    %wavelenght
    MPI_matrix = (gama_gap>1);
    TI_matrix = (gama_gap<1);
    delta_tilda_TI_gap=2/pi.*E_gap.*(sqrt(1+gama_gap.^2))./gama_gap.*ellipticE(1./(sqrt(1+gama_gap.^2)));%full keldish model in J from his paper 
    delta_tilda_mpi_gap=E_gap+(e^2.*F.^2)/(4.*m.*w.^2); %for the case of MPI from keldysh paper
    delta_tilda_gap = delta_tilda_TI_gap.*TI_matrix+delta_tilda_mpi_gap.*MPI_matrix; %[J]
    delta_tilda_ev_gap=delta_tilda_gap./e; %[eV]
    
    MPI_matrix = (gama_int>1);
    TI_matrix = (gama_int<1);
    delta_tilda_TI_int=2/pi.*E_int.*(sqrt(1+gama_int.^2))./gama_int.*ellipticE(1./(sqrt(1+gama_int.^2)));%full keldish model in J from his paper 
    delta_tilda_mpi_int=E_int+(e^2.*F.^2)/(4.*m.*w.^2); %for the case of MPI from keldysh paper
    delta_tilda_int = delta_tilda_TI_int.*TI_matrix+delta_tilda_mpi_int.*MPI_matrix; %[J]
    
    %% Keldysh strong feild ionization 
    %from Wavelength dependence of femtosecond laser-induced breakdown in water and implications Wavelength dependence of femtosecond laser-induced breakdown in water and implications
    dn_sfi_rho_init_c = n_sfi(delta_tilda_int,gama_int,w);

    [t_rho_init,rho_init]=ode45(@(t_solve,n)dn_sft_solver(t_solve,n,t,dn_sfi_rho_init_c,n_sat),t,0);

    dn_sfi_rho_gap = n_sfi(delta_tilda_gap,gama_gap,w);
    rho_init_transpose=transpose(rho_init);
    diff_n_sfi = dn_sfi_rho_init_c.*(1-rho_init_transpose/n_sat)+dn_sfi_rho_gap;

end
