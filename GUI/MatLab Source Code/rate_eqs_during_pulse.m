function dn_dt = rate_eqs_during_pulse(t, n, I_func, sigma_01, sigma_04, sigma_14, lifetime_1, lifetime_4, t_dephase_1, t_dephase_4, eps)
    try
        I = max(I_func(t), 0);  
        n0 = max(n(1), 0); n1 = max(n(2), 0); n1_d = max(n(3), 0); 
        n4 = max(n(4), 0); n4_d = max(n(5), 0); 
        lifetime_1 = max(lifetime_1, eps);
        lifetime_4 = max(lifetime_4, eps);
        t_dephase_1 = max(t_dephase_1, eps);
        t_dephase_4 = max(t_dephase_4, eps);
        dn0_dt = -I*(n0-n1)*sigma_01 - I^2/2*(n0-n4)*sigma_04;
        dn1_dt = I*(n0-n1)*sigma_01 - I*(n1-n4)*sigma_14 - n1/t_dephase_1 - n1/lifetime_1;
        dn1_d_dt = n1/t_dephase_1 - n1_d/lifetime_1 - I*(n1_d - n4)*sigma_14;
        dn4_dt = I^2/2*(n0-n4)*sigma_04 + I*(n1 - n4)*sigma_14 + I*(n1_d - n4)*sigma_14 - n4/t_dephase_4 - n4/lifetime_4;
        dn4_d_dt = n4/t_dephase_4 - n4_d/lifetime_4;
        dn01_dt = n1_d/lifetime_1 + n1/lifetime_1;
        dn04_dt = n4_d/lifetime_4 + n4/lifetime_4;
        dn_dt = [dn0_dt; dn1_dt; dn1_d_dt; dn4_dt; dn4_d_dt; dn01_dt; dn04_dt];
    catch
        dn_dt = zeros(7,1);  % Return zeros if any calculation fails 
    end
end