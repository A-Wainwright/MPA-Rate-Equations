% ==========================
% Dissociation ODE function
% ==========================
function dn_dt = dissociationODE(t, n, I_func, sigma_01, sigma_04, sigma_14, lifetime_1, lifetime_4, eps)
    try
        I = max(I_func(t), 0);
        n0 = n(1); n1 =n(2); n4 =n(3); n01 =n(4); n04=n(5); 

        lifetime_1 = max(lifetime_1, eps);
        lifetime_4 = max(lifetime_4, eps);

        dn0_dt = -I*(n0-n1)*sigma_01 - I^2/2*(n0-n4)*sigma_04;
        dn1_dt = I*(n0-n1)*sigma_01 - I*(n1-n4)*sigma_14 - n1/lifetime_1;
        dn4_dt = I^2/2*(n0-n4)*sigma_04 + I*(n1 - n4)*sigma_14 - n4/lifetime_4;
        dn01_dt = n1/lifetime_1;
        dn04_dt = n4/lifetime_4;
        dn_dt = [dn0_dt; dn1_dt; dn4_dt; dn01_dt; dn04_dt];
    catch
        dn_dt = zeros(5,1);  % Return zeros if any calculation fails 
    end
end

