function dn_dt = rate_eqs_after_pulse(~, n, lifetime_1, lifetime_4, t_dephase_1, t_dephase_4, eps)
    try
        n1 = max(n(2), 0); n1_d = max(n(3), 0); 
        n4 = max(n(4), 0); n4_d = max(n(5), 0); 
        % Set up nonzero denominators
        lifetime_1 = max(lifetime_1, eps);
        lifetime_4 = max(lifetime_4, eps);
        t_dephase_1 = max(t_dephase_1, eps);
        t_dephase_4 = max(t_dephase_4, eps);
        dn0_dt = 0;
        dn1_dt = -n1/t_dephase_1 - n1/lifetime_1;
        dn1_d_dt = n1/t_dephase_1 - n1_d/lifetime_1;
        dn4_dt = -n4/t_dephase_4 - n4/lifetime_4;
        dn4_d_dt = n4/t_dephase_4 - n4_d/lifetime_4;
        dn01_dt = n1_d/lifetime_1 + n1/lifetime_1;
        dn04_dt = n4_d/lifetime_4 + n4/lifetime_4;
        dn_dt = [dn0_dt; dn1_dt; dn1_d_dt; dn4_dt; dn4_d_dt; dn01_dt; dn04_dt];
    catch
        dn_dt = zeros(7,1);
    end
end