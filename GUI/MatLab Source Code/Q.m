
function [q] = Q(gamma, x)
    % Q function for Keldysh strong field ionization
    % Optimized version with reduced summation terms
    [K_gamma, E_gamma] = ellipke(gamma ./ (sqrt(1 + gamma.^2)));
    [K_val, E_val] = ellipke(1 ./ (sqrt(1 + gamma.^2)));
    
    % Reduced top value for faster computation (was 500, now adaptive)
    top = min(100, max(20, round(50./sqrt(gamma + 0.1)))); % Adaptive based on gamma
    x_val = 0:top;
    
    % Vectorized computation
    l_matrix = repmat(x_val(:), 1, length(x));
    x_matrix = repmat(x(:)', length(x_val), 1);
    
    exp_term = exp(-pi .* l_matrix .* (K_gamma - E_gamma) ./ E_val);
    dawson_arg = sqrt(pi .* (2 .* floor(x_matrix + 1) - 2 .* x_matrix + l_matrix) ./ (2 .* K_val .* E_val));
    dawson_term = dawson(dawson_arg);
    
    Q_sum_val = sum(exp_term .* dawson_term, 1);
    
    % Value being returned
    q = sqrt(pi ./ (2 .* K_val)) .* Q_sum_val;
end
