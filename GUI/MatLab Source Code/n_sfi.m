
function [n_gamma_factor_constant] = n_sfi(delta_tilda, gamma, w)
    % Strong field ionization rate calculation using Keldysh theory
    const = getPhysicalConstants();
    
    [K_1, E_1] = ellipke(gamma ./ (sqrt(1 + gamma.^2)));
    E_2 = ellipticE(1 ./ (sqrt(1 + gamma.^2))); % elliptic integral of the second kind
    x = delta_tilda ./ (const.h_bar .* w); % Number of photons needed to exceed the band gap
    Q_val = Q(gamma, x);

    % Full Keldysh model
    n_gamma_factor_constant = 2 .* w ./ (9 .* pi) .* ...
                              ((sqrt(1 + gamma.^2) ./ gamma) .* const.m .* w ./ const.h_bar).^(3/2) .* ...
                              Q_val .* exp(-pi .* floor(x + 1) .* (K_1 - E_1) ./ E_2);
end
