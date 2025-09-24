function checkCompilationReadiness()
    fprintf('Checking compilation readiness...\n');
    
    % Check for required toolboxes
    requiredToolboxes = {'MATLAB', 'Optimization_Toolbox'};
    installedToolboxes = ver;
    toolboxNames = {installedToolboxes.Name};
    
    for i = 1:length(requiredToolboxes)
        if any(contains(toolboxNames, requiredToolboxes{i}, 'IgnoreCase', true))
            fprintf('✓ %s found\n', requiredToolboxes{i});
        else
            fprintf('✗ %s NOT found - may cause compilation issues\n', requiredToolboxes{i});
        end
    end
    
    % Check file dependencies
    try
        [fList, pList] = matlab.codetools.requiredFilesAndProducts(mfilename('fullpath'));
        fprintf('\nRequired files (%d total):\n', length(fList));
        for i = 1:min(5, length(fList))  % Show first 5
            fprintf('  %s\n', fList{i});
        end
        if length(fList) > 5
            fprintf('  ... and %d more files\n', length(fList) - 5);
        end
        
        fprintf('\nRequired products (%d total):\n', length(pList));
        for i = 1:length(pList)
            fprintf('  %s\n', pList(i).Name);
        end
    catch ME
        fprintf('Warning: Could not analyze dependencies: %s\n', ME.message);
    end
    
    fprintf('\nCompilation readiness check complete!\n');
end