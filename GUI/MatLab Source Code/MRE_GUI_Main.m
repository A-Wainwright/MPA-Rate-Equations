function MRE_GUI_Main()
%% Filename: MPA_GUI_Main.m
% Author: Alexander Wainwright, Syeda Mahdia 
% Date Created: 2025-09-29
% Last Modified: 2025-09-29
% Version: 1.0 B
%
% Description:
%   GUI for modeling multi-photon absorption (MPA) in protein and water systems.
%   Allows users to input laser parameters, select molecular targets, and simulate
%   one-photon and two-photon absorption, strong-field ionization, and related dynamics.
%
% Inputs:
%   User inputs via GUI components:
%       - Laser wavelength (nm)
%       - Laser intensity (W/cm^2)
%       - Pulse duration (fs)
%       - Target molecules (proteins, water, metabolites)
%       - Solver options (ODE45, ODE23, ODE15s, Euler)
%       - Advanced parameters (cross-sections, excited state lifetimes)
%
% Outputs:
%   - Simulated population dynamics of molecular and electronic states
%   - Plots of absorption, ionization probability, and population over time
%   - Exportable data for further analysis
%
% Dependencies:
%   - MATLAB R2023a or later
%   - Signal Processing Toolbox
%   - Optimization Toolbox (if using advanced fitting)
%
% Usage:
%   Run `MPA_GUI_Main` in MATLAB to launch the GUI.
%
% Notes:
%   - Ensure laser pulse duration is compatible with excited state lifetimes.
%   - Warning alerts are included for intensities exceeding ionization thresholds.
%
% Revision History:
%   v1.0 - Initial beta testing release with GUI for MPA modeling, ODE solver options, and plotting.

%% 
if isdeployed
    fprintf('Running as compiled application\n');
else
    fprintf('Running in MATLAB environment\n');
end

% Layout constants
FIG_WIDTH = 860;
FIG_HEIGHT = 800;
TOP_PADDING = 25;
SIDE_PADDING = 20;
COL_SPACING = 40;
ROW_SPACING = 20;
LABEL_WIDTH = 250;
EDIT_WIDTH = 110;
FIELD_HEIGHT = 32;

% Color scheme
MAIN_BG = [0.98, 0.98, 0.98];
PANEL_BG = [0.95, 0.95, 0.95];
ACCENT_BG = [0.90, 0.90, 0.90];
BUTTON_BG = [0.25, 0.35, 0.45];
DARK_ACCENT = [0.15, 0.15, 0.15];
WHITE_TEXT = [1.0, 1.0, 1.0];

% Create main figure
fig = uifigure('Name','Multi-Photon Effects','Position',[100, 100, FIG_WIDTH, FIG_HEIGHT], ...
    'Color', MAIN_BG);

% Title
uilabel(fig,'Position',[20, FIG_HEIGHT-35, FIG_WIDTH-40, 30], ...
    'Text','MPA Excitation, Dissociation & Ionization', ...
    'FontSize',20,'FontWeight','bold','FontName','Segoe UI', ...
    'HorizontalAlignment','center','BackgroundColor',MAIN_BG,'FontColor',DARK_ACCENT);

% Process selection
uilabel(fig,'Position',[30, FIG_HEIGHT-70, 80, 28],'Text','PROCESS:', ...
    'FontWeight','bold','FontColor',DARK_ACCENT,'FontName','Segoe UI', ...
    'BackgroundColor',MAIN_BG,'FontSize',15);
processMenu = uidropdown(fig,'Items',{'Excitation','Dissociation','Ionization'}, ...
    'Position',[120, FIG_HEIGHT-70, 150, 28],'Value','Excitation', ...
    'BackgroundColor',[1,1,1],'FontColor',DARK_ACCENT,'FontName','Segoe UI','FontSize',13);

% Sample name input
uilabel(fig,'Position',[300, FIG_HEIGHT-70, 140, 28],'Text','Name of Sample:', ...
    'FontWeight','bold','FontColor',DARK_ACCENT,'FontName','Segoe UI', ...
    'BackgroundColor',MAIN_BG,'FontSize',15);
uieditfield(fig,'text','Tag','sampleName','Value','My Sample', ...
    'Position',[450, FIG_HEIGHT-70, 180, 28],'BackgroundColor',[1,1,1], ...
    'FontColor',DARK_ACCENT,'FontName','Segoe UI','FontSize',13);

% Parameter panel
paramPanel = uipanel(fig,'Title','Input Variables','Position',[20,170,FIG_WIDTH-40,FIG_HEIGHT-250], ...
    'BackgroundColor',PANEL_BG,'ForegroundColor',DARK_ACCENT,'FontWeight','bold', ...
    'FontName','Segoe UI','FontSize',13);

% Solver Panel
solverGroup = uibuttongroup(fig,'Title','ODE Solver','Position',[20,20,FIG_WIDTH-40,100], ...
    'BackgroundColor',PANEL_BG,'FontName','Segoe UI','FontSize',13,'FontWeight','bold');

uiradiobutton(solverGroup,'Text','ode45 (non-stiff)','Position',[10,55,200,20], ...
    'FontName','Segoe UI','FontSize',12,'Tag','ode45'); % default

uiradiobutton(solverGroup,'Text','ode23 (non-stiff)','Position',[160,55,200,20], ...
    'FontName','Segoe UI','FontSize',12,'Tag','ode23');

uiradiobutton(solverGroup,'Text','ode15s (stiff)','Position',[310,55,200,20], ...
    'FontName','Segoe UI','FontSize',12,'Tag','ode15s');

uiradiobutton(solverGroup,'Text','ode23s (stiff)','Position',[460,55,200,20], ...
    'FontName','Segoe UI','FontSize',12,'Tag','ode23s');

uiradiobutton(solverGroup,'Text','Euler (explicit)','Position',[610,55,200,20], ...
    'FontName','Segoe UI','FontSize',12,'Tag','Euler');

% Run button
uibutton(fig,'Text','Run Simulation','Position',[(FIG_WIDTH-140)/2-200,30,140,40], ...
    'ButtonPushedFcn',@(btn,event) onRun(processMenu.Value), ...
    'BackgroundColor',BUTTON_BG,'FontColor',WHITE_TEXT,'FontWeight','bold', ...
    'FontName','Segoe UI','FontSize',14);

% Close all plots
uibutton(fig, 'Text', 'Close Plots', 'Position', [(FIG_WIDTH-140)/2+200, 30, 140, 40], ...
    'ButtonPushedFcn', @(btn,event) Close_all(), ...
    'BackgroundColor', BUTTON_BG, ...
    'FontColor', WHITE_TEXT, ...
    'FontWeight', 'bold', ...
    'FontName', 'Segoe UI', ...
    'FontSize', 14);

    function Close_all()
        close all
    end

% Set up callbacks and initialize
processMenu.ValueChangedFcn = @(dd,event) updateInputs(dd.Value);
advancedPanel = [];
updateInputs('Excitation');

%% --- Nested function: updateInputs ---
    function updateInputs(mode)
        % Clear existing components in paramPanel
        delete(paramPanel.Children);

        switch mode
            case 'Excitation'
                createInput(paramPanel,'Protein Concentration (M)','conc',30e-3);
                createInput(paramPanel,'Absorption Coefficient α (cm^{-1})','alpha',1e4);
                createInput(paramPanel,'σ_{01} (cm^2)','sigma',1e-20);
                createInput(paramPanel,'Photons per Chromophore (Calc.)','phot_per_chrom',0);

                createInput(paramPanel,'Wavelength (nm)','wavelength',530);
                createInput(paramPanel,'Pulse Duration (s)','tau_p',200e-15);
                createInput(paramPanel,'Laser Peak Power (W/cm^2)','I',1e13);

                createInput(paramPanel,'σ_{0n} (m^4·s/photon)','sigma02',290e-58);
                createInput(paramPanel,'σ_{1n} (cm^2)','sigma12',1e-20);
                createInput(paramPanel,'Dephasing Time S1 (s)','t_dephase_1',28e-15);
                createInput(paramPanel,'Decay Lifetime S1 (s)','lifetime_1',450e-15);
                createInput(paramPanel,'Decay Lifetime S4 (s)','lifetime_4',200e-15);
                createInput(paramPanel,'Max depth (m)','z_max',3e-6);

                solverGroup.SelectedObject = findobj(solverGroup,'Tag','ode45'); % default: ode45

            case 'Dissociation'
                createInput(paramPanel,'Protein Concentration (M)','conc',30e-3);
                createInput(paramPanel,'Absorption Coefficient α (cm^{-1})','alpha',1e4);
                createInput(paramPanel,'σ_{01} (cm^2)','sigma',1e-20);
                createInput(paramPanel,'Photons per Chromophore (Calc.)','phot_per_chrom',0);

                createInput(paramPanel,'σ_{0n} (m^4·s/photon)','sigma02',1e-50);
                createInput(paramPanel,'σ_{1n} (m^2)','sigma12',1e-20);
                createInput(paramPanel,'τ_p (s)','tau_p',200*10^-15);
                createInput(paramPanel,'Peak Intensity (W/cm^2)','I',1e12);
                createInput(paramPanel,'Sample Depth (μm)','d_sample',5);
                createInput(paramPanel,'Relaxation Time τ_1 (s)','tau_1',1000*10^-15);
                createInput(paramPanel,'Relaxation Time τ_2 (s)','tau_2',800*10^-15);
                createInput(paramPanel,'Wavelength (nm)','wavelength',800);

                solverGroup.SelectedObject = findobj(solverGroup,'Tag','Euler'); % default: ode15s

                createAdvancedDissociationInputs();

            case 'Ionization'
                createInput(paramPanel,'Wavelength (m)','wavelength',515e-9);
                createInput(paramPanel,'Pulse Duration (s)','tau_p',250e-15);
                createInput(paramPanel,'Peak Power FWHM (W/cm^2)','Peak_power_FWHM',6e12);
                createInput(paramPanel,'Bandgap Energy (eV)','E_gap',9.5);
                createInput(paramPanel,'Intermediate Level Energy (eV)','E_int',6.6);
                createInput(paramPanel,'Electron Collision Time (s)','t_coll',0.9e-15);
                createInput(paramPanel,'Bound Electron Density (1/m^3)','n_bound',6.68e28);
                createInput(paramPanel,'Initial Neutral Density (cm^{-3})','n0_cm3',1e19);

                solverGroup.SelectedObject = findobj(solverGroup,'Tag','ode45'); % default: Euler
        end
    end

%% --- updateAbsorptionFields (with photons per chromophore) ---
    function updateAbsorptionFields(fig, src)
        NA = 6.022e23; % Avogadro number
        h = 6.626e-34; % Planck constant
        c = 3e8;       % speed of light

        % Find relevant fields
        concField  = findobj(fig,'Tag','conc');
        alphaField = findobj(fig,'Tag','alpha');
        sigmaField = findobj(fig,'Tag','sigma');
        photField  = findobj(fig,'Tag','phot_per_chrom');
        wavelengthField = findobj(fig,'Tag','wavelength'); % in nm

        % Labels for notes
        concLabel  = findobj(fig,'Tag','conc_note');
        alphaLabel = findobj(fig,'Tag','alpha_note');
        sigmaLabel = findobj(fig,'Tag','sigma_note');
        photLabel  = findobj(fig,'Tag','phot_note');

        if isempty(concField) || isempty(alphaField) || isempty(sigmaField) || isempty(photField)
            return;
        end

        C = concField.Value;      % M
        alpha = alphaField.Value; % cm^-1
        sigma = sigmaField.Value; % cm^2
        lambda_nm = wavelengthField.Value; % nm

        % Number density [cm^-3]
        N = C * NA * 1e-3;

        % Update sigma or alpha depending on which changed
        switch src.Tag
            case {'conc','alpha'}
                if N > 0
                    sigmaField.Value = str2double(sprintf('%.4g', alpha / N));
                    sigmaLabel.Text = 'σ = α / N';
                    alphaLabel.Text = 'User Value';
                else
                    sigmaField.Value = NaN;
                end
            case 'sigma'
                alphaField.Value = str2double(sprintf('%.4g', sigma * N));
                alphaLabel.Text = 'α = σ · N';
                sigmaLabel.Text = 'User Value';
        end

        % Update Photons per Chromophore
        if ~isempty(photField) && ~isempty(wavelengthField) && lambda_nm > 0
            I_field = findobj(fig,'Tag','I');       % Laser intensity (W/cm^2)
            tau_field = findobj(fig,'Tag','tau_p'); % Pulse duration (s)
            if ~isempty(I_field) && ~isempty(tau_field)
                I = I_field.Value;        % W/cm^2
                tau = tau_field.Value;    % s
                lambda_m = lambda_nm * 1e-9;
                photonEnergy = h*c/lambda_m;       % J per photon
                photonFlux = (I * tau) / photonEnergy; % photons/cm^2 per pulse
                sigma_cm2 = sigmaField.Value;          % cm^2
                photonsPerChrom = sigma_cm2 * photonFlux;

                photField.Value = str2double(sprintf('%.4g', photonsPerChrom));
                if ~isempty(photLabel)
                    photLabel.Text = 'Updated automatically';
                end

                % 🚨 LIVE WARNING CHECK
                if photonsPerChrom > 1
                    uialert(fig, ...
                        sprintf('Warning: Photons per chromophore = %.3g (greater than 1).', photonsPerChrom), ...
                        'Photon Warning', ...
                        'Icon','warning');
                end

                if strcmp(src.Tag,'I') && src.Value > 1e13
                    uialert(fig, ...
                        sprintf(['Warning: An intensity of %.2e W/cm^2 risks ionization.\n' ...
                        'Please run the ionization model to assess the risk of ionization with your input parameters.'], src.Value), ...
                        'High Intensity Warning', ...
                        'Icon','warning');
                end
            end
        end

        % Format fields
        concField.Value  = str2double(sprintf('%.4g', concField.Value));
        alphaField.Value = str2double(sprintf('%.4g', alphaField.Value));
        sigmaField.Value = str2double(sprintf('%.4g', sigmaField.Value));
        photField.Value  = str2double(sprintf('%.4g', photField.Value));
    end


%% --- createInput ---
    function createInput(parent,labelText,tag,defaultVal,isAdvanced)
        if nargin<5, isAdvanced=false; end

        LABEL_WIDTH = 250; EDIT_WIDTH = 110; FIELD_HEIGHT=32; SIDE_PADDING=20; ROW_SPACING=10;
        ACCENT_BG = [0.90 0.90 0.90]; PANEL_BG=parent.BackgroundColor; DARK_ACCENT=[0.15 0.15 0.15];

        numExistingItems = sum(arrayfun(@(x) isa(x,'matlab.ui.control.Label'), parent.Children));
        rowIdx = floor(numExistingItems/2); colIdx = mod(numExistingItems,2);
        panelHeight = parent.Position(4);
        yPos = panelHeight - (rowIdx+1)*FIELD_HEIGHT - rowIdx*ROW_SPACING-50;
        xLabel = SIDE_PADDING + colIdx*(LABEL_WIDTH+EDIT_WIDTH+SIDE_PADDING);
        xEdit  = xLabel + LABEL_WIDTH;

        % Label
        uilabel(parent,'Text',labelText,'Position',[xLabel,yPos,LABEL_WIDTH,FIELD_HEIGHT], ...
            'FontSize',isAdvanced*19+(~isAdvanced)*15,'FontWeight','bold','FontName','Segoe UI', ...
            'FontColor',DARK_ACCENT,'BackgroundColor',isAdvanced*ACCENT_BG+(~isAdvanced)*PANEL_BG, ...
            'Interpreter','tex');

        % Numeric edit field
        editField = uieditfield(parent,'numeric','Tag',tag,'Value',defaultVal, ...
            'Position',[xEdit,yPos,EDIT_WIDTH,FIELD_HEIGHT],'FontSize',isAdvanced*19+(~isAdvanced)*15, ...
            'FontName','Segoe UI','FontColor',DARK_ACCENT,'BackgroundColor',[1 1 1]);

        % Note label for conc/alpha/sigma/phot
        if ismember(tag,{'conc','alpha','sigma','phot_per_chrom'})
            uilabel(parent,'Tag',[tag '_note'],'Text','','Position',[xEdit+EDIT_WIDTH+10,yPos,150,FIELD_HEIGHT], ...
                'FontSize',12,'FontColor',[0.2 0.2 0.2],'BackgroundColor',parent.BackgroundColor);
        end

        % ValueChangedFcn
        editField.ValueChangedFcn = @(src,event) formatAndUpdate(src,parent);

        function formatAndUpdate(src,parentPanel)
            src.Value = str2double(sprintf('%.4g',src.Value));
            if ismember(src.Tag,{'conc','alpha','sigma','phot_per_chrom','I','tau_p','wavelength'})
                try
                    updateAbsorptionFields(parentPanel.Parent,src);
                catch
                end
            end

            % Check for pulse duration exceeding excited-state lifetime
            tauField = findobj(fig,'Tag','tau_p');       % Pulse duration (s)
            lifetimeField = findobj(fig,'Tag','lifetime_1'); % Excited-state lifetime (s)

            if ~isempty(tauField) && ~isempty(lifetimeField)
                tau = tauField.Value;
                lifetime = lifetimeField.Value;

                if tau > lifetime
                    uialert(fig, ...
                        sprintf(['Warning: Pulse duration (%.2e s) exceeds the excited-state lifetime (%.2e s).\n' ...
                        'This model is not suitable for such a system. Please reduce the pulse duration.'], tau, lifetime), ...
                        'Pulse Duration Warning', ...
                        'Icon','warning');
                end
            end

        end
    end

%% --- createAdvancedDissociationInputs ---
    function createAdvancedDissociationInputs()
        % Creates the extra panel for dissociation
        num_main_inputs = numel(paramPanel.Children) / 2;
        panelHeight = paramPanel.Position(4);
        advanced_y_pos = panelHeight - TOP_PADDING - (ceil(num_main_inputs/2) + 1) * FIELD_HEIGHT - ceil(num_main_inputs/2) * ROW_SPACING-70;

        advancedToggle = uibutton(paramPanel, 'push', ...
            'Text', '▼ Advanced Inputs', ...
            'Position', [SIDE_PADDING, advanced_y_pos, 180, 30], ...
            'ButtonPushedFcn', @toggleAdvancedInputs, ...
            'BackgroundColor', BUTTON_BG, ...
            'FontColor', WHITE_TEXT, ...
            'FontWeight', 'bold', ...
            'FontName', 'Segoe UI', ...
            'FontSize', 17);

        advancedPanel = uipanel(paramPanel, ...
            'Title', '', ...
            'Position', [SIDE_PADDING, advanced_y_pos - 100, FIG_WIDTH - 80, 90], ...
            'BorderType', 'line', ...
            'BackgroundColor', ACCENT_BG);

        createInput(advancedPanel, 'Spatial Resolution N_z', 'Nz', 100, true);

        % Store references
        advancedPanel.UserData = struct('toggleButton', advancedToggle, 'isExpanded', true);

        % 🔑 Bring the panel on top of all siblings
        uistack(advancedPanel, 'top');
    end


%% --- toggleAdvancedInputs ---
    function toggleAdvancedInputs(src,~)
        targetPanel = [];
        parentChildren = src.Parent.Children;
        for i=1:length(parentChildren)
            if isa(parentChildren(i),'matlab.ui.container.Panel') && ...
                    isfield(parentChildren(i).UserData,'toggleButton') && ...
                    parentChildren(i).UserData.toggleButton == src
                targetPanel = parentChildren(i);
                break;
            end
        end
        if isempty(targetPanel)||~isvalid(targetPanel), return; end
        toggleData = targetPanel.UserData;
        if toggleData.isExpanded
            targetPanel.Visible = 'off';
            src.Text = strrep(src.Text,'▼','▶');
            toggleData.isExpanded = false;
        else
            targetPanel.Visible = 'on';
            src.Text = strrep(src.Text,'▶','▼');
            toggleData.isExpanded = true;
        end
        targetPanel.UserData = toggleData;
    end

%% --- onRun ---
    function onRun(mode)
        try
            inputs = struct();
            sampleNameField = findobj(fig,'Tag','sampleName');
            inputs.sampleName = sampleNameField.Value;
            if isempty(inputs.sampleName)
                inputs.sampleName = 'DefaultSample';
            end

            selectedSolver = solverGroup.SelectedObject.Tag;
            inputs.solver = selectedSolver;

            for field = paramPanel.Children'
                if isa(field,'matlab.ui.control.NumericEditField')
                    val = field.Value;
                    if isnan(val)||~isfinite(val)
                        uialert(fig,sprintf('Invalid value for %s',field.Tag),'Input Error');
                        return;
                    end
                    inputs.(field.Tag) = val;
                end
            end

            if ~isempty(advancedPanel) && isvalid(advancedPanel)
                for field = advancedPanel.Children'
                    if isa(field,'matlab.ui.control.NumericEditField')
                        val = field.Value;
                        if isnan(val)||~isfinite(val)
                            uialert(fig,sprintf('Invalid value for %s',field.Tag),'Input Error');
                            return;
                        end
                        inputs.(field.Tag) = val;
                    end
                end
            end

            % 🚨CHECK: photons per chromophore > 1
            if isfield(inputs,'phot_per_chrom') && inputs.phot_per_chrom > 1
                uialert(fig, ...
                    sprintf('Warning: Photons per chromophore = %.3g (greater than 1).', inputs.phot_per_chrom), ...
                    'Photon Warning', ...
                    'Icon','warning');
            end
            % 🚨CHECK Ionization: I > 1e13
            if isfield(inputs,'I')
                checkPowerWarning(mode, inputs.I);
            end

            % Run simulation (placeholders)
            switch mode
                case 'Excitation'
                    runExcitation(inputs);
                case 'Dissociation'
                    runDissociation(inputs);
                case 'Ionization'
                    runIonization(inputs);
            end

        catch ME
            uialert(fig,['Unexpected error: ',ME.message],'Error');
            fprintf('Error details: %s\n',getReport(ME));
        end
    end

    function checkPowerWarning(simType, powerValue)
        % Threshold in W/cm^2
        threshold = 1e13;

        if (strcmp(simType, 'Excitation') || strcmp(simType, 'Dissociation')) ...
                && powerValue > threshold
            warndlg(sprintf(['Warning: Selected power (%.2e W/cm^2) is high enough to cause signficant ioniation. ' ...
                'We recommended you run the ionization model as well as the %s model for I>1e13 W/cm^2.'], ...
                powerValue, simType), ...
                'High Power Warning');
        end
    end
end