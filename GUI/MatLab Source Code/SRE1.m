
function dn_dt = SRE1(t_solve, n, t, n_sfi_1, n_AI_Asym_1, n_bound)
    % Rate equation for total electron density evolution
    
    % Use faster interpolation with bounds checking
    if t_solve < t(1)
        n_sfi = n_sfi_1(1);
        n_AI_Asym = n_AI_Asym_1(1);
    elseif t_solve > t(end)
        n_sfi = n_sfi_1(end);
        n_AI_Asym = n_AI_Asym_1(end);
    else
        n_sfi = interp1(t, n_sfi_1, t_solve, 'linear');   
        n_AI_Asym = interp1(t, n_AI_Asym_1, t_solve, 'linear'); 
    end
    
    dn_dt = n_sfi + (n_AI_Asym) .* n .* (n_bound - n) ./ n_bound - 1.8e-9 * 1e-6 .* n.^2; % added recombination 
end
