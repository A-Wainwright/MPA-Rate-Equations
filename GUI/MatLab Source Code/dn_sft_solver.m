
function dndt = dn_sft_solver(t, n, t_sln, n_sfi_E_ini, rho_ini_max)
    % Rate equation solver for strong field ionization
    n_sfi_E_ini_val = interp1(t_sln, n_sfi_E_ini, t, 'linear', 0);   
    dndt = n_sfi_E_ini_val .* (1 - n / rho_ini_max);
end
