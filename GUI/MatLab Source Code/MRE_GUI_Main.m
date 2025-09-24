function MRE_GUI()
    
    % Check if running as compiled application
    if isdeployed
        % Running as compiled application - no need to modify path
        fprintf('Running as compiled application\n');
    else
        % Running in MATLAB environment
        fprintf('Running in MATLAB environment\n');
    end
    
    % Check compilation readiness
    checkCompilationReadiness();

    % Layout constants
    FIG_WIDTH = 860;  
    FIG_HEIGHT = 600; 
    TOP_PADDING = 25;
    SIDE_PADDING = 20;
    COL_SPACING = 40;  
    ROW_SPACING = 20;  
    LABEL_WIDTH = 260; 
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
    fig = uifigure('Name','Multi-Photon Effects Simulator','Position',[100, 100, FIG_WIDTH, FIG_HEIGHT], ...
                   'Color', MAIN_BG);
    
    % Title
    titleLabel = uilabel(fig,'Position',[20, FIG_HEIGHT-35, FIG_WIDTH-40, 30], ... 
                        'Text','MPA Excitation, Dissociation & Ionization Simulator', ...
                        'FontSize', 20, 'FontWeight', 'bold', ...
                        'FontName', 'Segoe UI', ...
                        'HorizontalAlignment', 'center', ...
                        'BackgroundColor', MAIN_BG, ...
                        'FontColor', DARK_ACCENT);
    
    % Process selection
    uilabel(fig,'Position',[30, FIG_HEIGHT-70, 80, 28],'Text','PROCESS:', ... 
           'FontWeight', 'bold', 'FontColor', DARK_ACCENT, ...
           'FontName', 'Segoe UI', ...
           'BackgroundColor', MAIN_BG, 'FontSize', 15); 
    processMenu = uidropdown(fig, 'Items',{'Excitation','Dissociation','Ionization'}, ...
                           'Position', [120, FIG_HEIGHT-70, 150, 28], 'Value', 'Excitation', ... 
                           'BackgroundColor', [1, 1, 1], 'FontColor', DARK_ACCENT, ...
                           'FontName', 'Segoe UI', ...
                           'FontSize', 13); 
    
    % Sample name input
    uilabel(fig, 'Position', [300, FIG_HEIGHT-70, 140, 28], 'Text', 'Name of Sample:', ...
           'FontWeight', 'bold', 'FontColor', DARK_ACCENT, ...
           'FontName', 'Segoe UI', 'BackgroundColor', MAIN_BG, 'FontSize', 15);
    uieditfield(fig, 'text', 'Tag', 'sampleName', 'Value', 'My Sample', ...
               'Position', [450, FIG_HEIGHT-70, 180, 28], ...
               'BackgroundColor', [1, 1, 1], 'FontColor', DARK_ACCENT, ...
               'FontName', 'Segoe UI', 'FontSize', 13);
    
    % Parameter panel
    paramPanel = uipanel(fig, 'Title', 'Input Variables', ...
                       'Position', [20, 85, FIG_WIDTH - 40, FIG_HEIGHT - 155], ... 
                       'BackgroundColor', PANEL_BG, ...
                       'ForegroundColor', DARK_ACCENT, ...
                       'FontWeight', 'bold', ...
                       'FontName', 'Segoe UI', ...
                       'FontSize', 13);
    
    % Run button
    uibutton(fig, 'Text', 'Run Simulation', 'Position', [(FIG_WIDTH-140)/2-200, 30, 140, 40], ...
             'ButtonPushedFcn', @(btn,event) onRun(processMenu.Value), ...
             'BackgroundColor', BUTTON_BG, ...
             'FontColor', WHITE_TEXT, ...
             'FontWeight', 'bold', ...
             'FontName', 'Segoe UI', ...
             'FontSize', 14); 
    
    % Close All Plots
    uibutton(fig, 'Text', 'Close Plots', 'Position', [(FIG_WIDTH-140)/2+200, 30, 140, 40], ...
             'ButtonPushedFcn', @(btn,event) Close_all(), ...
             'BackgroundColor', BUTTON_BG, ...
             'FontColor', WHITE_TEXT, ...
             'FontWeight', 'bold', ...
             'FontName', 'Segoe UI', ...
             'FontSize', 14); 
    
    % Set up callbacks and initialize
    processMenu.ValueChangedFcn = @(dd,event) updateInputs(dd.Value);
    advancedPanel = [];
    updateInputs('Excitation');
    function Close_all()
        close all
    end 
    % Nested function: updateInputs
    function updateInputs(mode)
        % Clear existing components in paramPanel
        delete(paramPanel.Children);
        switch mode
            case 'Excitation'
                createInput(paramPanel,'Wavelength (nm)', 'wavelength',530);
                createInput(paramPanel,'Laser Peak Power (W/cm^2)', 'power', 1e13);
                createInput(paramPanel,'Protein Molarity (M)','Mprot',30e-3);
                createInput(paramPanel, 'Extinction (1/(M·cm))', 'molar_extinction_Protein',49000);
                createInput(paramPanel,'CrossSection0_4 (m^4·s/photon)', 'Cross_section_0_4_per_photon', 290e-58); 
                createInput(paramPanel, 'Cross_section_1_4 (m^2/photon)','Cross_section_1_4_per_photon', 0);
                createInput(paramPanel, 'Dephasing Time S1 (s)','t_dephase_1', 28e-15);
                createInput(paramPanel, 'Decay Lifetime S1 (s)', 'lifetime_1', 450e-15);
                createInput(paramPanel, 'Decay Lifetime S4 (s)', 'lifetime_4', 200e-15);
                createInput(paramPanel, 'Max depth (m)','z_max',3e-6);
                createInput(paramPanel, 'Pulse Duration (FWHM,s)', 'tau_p', 200e-15);
            case 'Dissociation'
                createInput(paramPanel, 'Cross Section σ_{01} (m^2)', 'sigma01', 1e-20);
                createInput(paramPanel, 'Cross Section σ_{02} (m^4·s/photon)', 'sigma02', 1e-50);  
                createInput(paramPanel, 'Cross Section σ_{12} (m^2)', 'sigma12', 1e-20);
                createInput(paramPanel, 'Pulse Duration τ_p (fs)', 'tau_p', 200);
                createInput(paramPanel, 'Peak Intensity (W/cm^2)', 'I', 1e12); 
                createInput(paramPanel, 'Sample Depth (μm)', 'd_sample', 5);
                createInput(paramPanel, 'Relaxation Time τ_1 (fs)', 'tau_1', 1000); 
                createInput(paramPanel, 'Relaxation Time τ_2 (fs)', 'tau_2', 800);
                createInput(paramPanel, 'Wavelength (nm)', 'wavelength', 800);
                createInput(paramPanel, 'Absorption Coefficient α (m^{-1})', 'alpha', 1e4);
                
                createAdvancedDissociationInputs();
               
            case 'Ionization' 
                createInput(paramPanel, 'Wavelength (m)', 'wavelength', 515e-9); 
                createInput(paramPanel, 'Pulse Duration (s)', 'tau_p', 250e-15); 
                createInput(paramPanel, 'Spot Size (μm)', 'omega_um', 2.5); 
                createInput(paramPanel, 'Peak Power FWHM (W/cm²)', 'Peak_power_FWHM', 6e12); 
                createInput(paramPanel, 'Bandgap Energy (eV)', 'E_gap', 9.5); 
                createInput(paramPanel, 'Intermediate Level Energy (eV)', 'E_int', 6.6);
                createInput(paramPanel, 'Electron Collision Time (s)', 't_coll', 0.9e-15); 
                createInput(paramPanel, 'Bound Electron Density (1/m³)', 'n_bound', 6.68e28); 
                createInput(paramPanel, 'Initial Neutral Density (cm⁻³)', 'n0_cm3', 1e19); 
        end
    end

    % Nested function: createAdvancedDissociationInputs
    function createAdvancedDissociationInputs()
        % Creates the extra panel for dissociation
        num_main_inputs = numel(paramPanel.Children) / 2;
        panelHeight = paramPanel.Position(4);
        advanced_y_pos = panelHeight - TOP_PADDING - (ceil(num_main_inputs/2) + 1) * FIELD_HEIGHT - ceil(num_main_inputs/2) * ROW_SPACING - 50;
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
            'Position', [SIDE_PADDING, advanced_y_pos - 70, FIG_WIDTH - 80, 65], ...
            'BorderType', 'line', ...
            'BackgroundColor', ACCENT_BG);
        createInput(advancedPanel, 'Spatial Resolution N_z', 'Nz', 100, true); 
        advancedPanel.UserData = struct('toggleButton', advancedToggle, 'isExpanded', true);
    end

    % Nested function: toggleAdvancedInputs
    function toggleAdvancedInputs(src, ~)
        if isfield(src.UserData, 'associatedPanel')
            targetPanel = src.UserData.associatedPanel;
        else
            targetPanel = [];
            parentChildren = src.Parent.Children;
            for i = 1:length(parentChildren)
                if isa(parentChildren(i), 'matlab.ui.container.Panel') && ...
                   isfield(parentChildren(i).UserData, 'toggleButton') && ...
                   parentChildren(i).UserData.toggleButton == src
                    targetPanel = parentChildren(i);
                    break;
                end
            end
        end
        if isempty(targetPanel) || ~isvalid(targetPanel)
            return;
        end
        toggleData = targetPanel.UserData;
        if toggleData.isExpanded
            % Collapse
            targetPanel.Visible = 'off';
            src.Text = strrep(src.Text, '▼', '▶');
            toggleData.isExpanded = false;
        else
            % Expand
            targetPanel.Visible = 'on';
            src.Text = strrep(src.Text, '▶', '▼');
            toggleData.isExpanded = true;
        end
        targetPanel.UserData = toggleData;
    end

    % Nested function: createInput
    function createInput(parent, label, tag, defaultVal, isAdvanced)
        if nargin < 5
            isAdvanced = false;
        end
        try
            if isAdvanced
                uilabel(parent, 'Text', label, 'Position', [10, 20, LABEL_WIDTH, FIELD_HEIGHT], ... 
                       'FontColor', DARK_ACCENT, 'BackgroundColor', ACCENT_BG, ...
                       'FontWeight', 'bold', 'FontName', 'Segoe UI', ...
                       'FontSize', 19);
                uieditfield(parent, 'numeric', 'Tag', tag, 'Value', defaultVal, ...
                            'Position', [LABEL_WIDTH + 15, 20, EDIT_WIDTH, FIELD_HEIGHT], ...
                            'BackgroundColor', [1, 1, 1], 'FontColor', DARK_ACCENT, ...
                            'FontName', 'Segoe UI', ...
                            'FontSize', 19); 
            else
                num_existing_items = numel(parent.Children) / 2;
                row_idx = floor(num_existing_items / 2);
                col_idx = mod(num_existing_items, 2);
                panelHeight = parent.Position(4);
                y_pos = panelHeight - TOP_PADDING - (row_idx + 1) * FIELD_HEIGHT - row_idx * ROW_SPACING;
                x_label_pos = SIDE_PADDING + col_idx * (LABEL_WIDTH + EDIT_WIDTH + COL_SPACING);
                x_edit_pos = x_label_pos + LABEL_WIDTH;
                uilabel(parent, 'Text', label, 'Position', [x_label_pos, y_pos, LABEL_WIDTH, FIELD_HEIGHT], ...
                       'FontColor', DARK_ACCENT, 'BackgroundColor', PANEL_BG, ...
                       'FontWeight', 'bold', 'FontName', 'Segoe UI', ...
                       'FontSize', 15); 
                uieditfield(parent, 'numeric', 'Tag', tag, 'Value', defaultVal, ...
                            'Position', [x_edit_pos, y_pos, EDIT_WIDTH, FIELD_HEIGHT], ...
                            'BackgroundColor', [1, 1, 1], 'FontColor', DARK_ACCENT, ...
                            'FontName', 'Segoe UI', ...
                            'FontSize', 15); 
            end
        catch ME
            warning('Failed to create input field %s: %s', tag, ME.message);
        end
    end

    % Nested function: onRun
    function onRun(mode)
        try
            % Input parsing
            inputs = struct();
            sampleNameField = findobj(fig, 'Tag', 'sampleName');
            inputs.sampleName = sampleNameField.Value;
            if isempty(inputs.sampleName)
                inputs.sampleName = 'DefaultSample';
            end
            
            for field = paramPanel.Children'
                if isa(field, 'matlab.ui.control.NumericEditField')
                    val = field.Value;
                    if isnan(val) || ~isfinite(val)
                        uialert(fig, sprintf('Invalid value for %s. Please enter a valid number.', field.Tag), 'Input Error');
                        return;
                    end
                    inputs.(field.Tag) = val;
                end
            end
           
            if strcmp(mode, 'Dissociation') && ~isempty(advancedPanel) && isvalid(advancedPanel)
                for field = advancedPanel.Children'
                    if isa(field, 'matlab.ui.control.NumericEditField')
                        val = field.Value;
                        if isnan(val) || ~isfinite(val)
                            uialert(fig, sprintf('Invalid value for %s. Please enter a valid number.', field.Tag), 'Input Error');
                            return;
                        end
                        inputs.(field.Tag) = val;
                    end
                end
            end
            
            if ~validateInputs(inputs, mode)
                return;
            end
            
            % Run simulation with error handling
            switch mode
                case 'Excitation'
                    try
                        runExcitation(inputs);
                    catch ME
                        uialert(fig, ['Error in Excitation simulation: ' ME.message], 'Simulation Error');
                        fprintf('Excitation error details: %s\n', getReport(ME));
                    end
                case 'Dissociation'
                    try
                        runDissociation(inputs);
                    catch ME
                        uialert(fig, ['Error in Dissociation simulation: ' ME.message], 'Simulation Error');
                        fprintf('Dissociation error details: %s\n', getReport(ME));
                    end
                case 'Ionization'
                     try
                        runIonization(inputs);
                    catch ME
                      uialert(fig, ['Error in Ionization simulation: ' ME.message], 'Simulation Error');
                      fprintf('Ionization error details: %s\n', getReport(ME));
                    end
            end
        catch ME
            uialert(fig, ['Unexpected error: ' ME.message], 'Error');
            fprintf('Unexpected error details: %s\n', getReport(ME));
        end
    end

    % Nested function: validateInputs
    function isValid = validateInputs(inputs, mode)
        isValid = true;
        % Common validations
        fieldNames = fieldnames(inputs);
        for i = 1:length(fieldNames)
            % Skip validation for sampleName
            if strcmp(fieldNames{i}, 'sampleName'), continue; end
            val = inputs.(fieldNames{i});
            if val <= 0 && ~strcmp(fieldNames{i}, 'Cross_section_1_4_per_photon')
                uialert(fig, sprintf('%s must be positive and greater than 0', fieldNames{i}), 'Input Error');
                isValid = false;
                return;
            end
        end
        % Mode-specific validation
        switch mode
            case 'Excitation'
                if inputs.wavelength < 100 || inputs.wavelength > 2000
                    uialert(fig, 'Wavelength should be between 100-2000 nm', 'Input Error');
                    isValid = false;
                end
                if inputs.power > 1e16
                    uialert(fig, 'Power too high - may cause numerical instability', 'Input Warning');
                end
   
            case 'Dissociation'
                if inputs.tau_p > 10000
                    uialert(fig, 'Pulse duration too large (>10 ps)', 'Input Warning');
                end
        end
    end

end
