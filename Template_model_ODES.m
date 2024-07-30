function dydt = Template_model_ODES(t, y, p)
    % Parameters
    DEP = 1.3;
    DEN = 0.7;
    vBS = p(1) * DEN;
    v_pl1 = p(2) * DEP;
    v_pl3 = p(3);
    v_pt = p(4) * DEP;
    v_ptt = p(5) * DEN;
    v_es = p(6); 
    v_d = p(7) * DEP;
    v_s = p(8) * DEP;
    v_p = p(9)* DEN;
    v_pl2 = p(10)* DEP;
    v_p14 = p(11) * DEN; 
    DietL = p(12); 
    a3 = p(13);
    a4 = p(14);
   
    
    % Extract state variables from y
    Chep = y(1);
    Cstor = y(2);
    Cpt = y(3);
    Cpla = y(4);
    Ces = y(5);


    % Define the system of ODEs
    dxdt = zeros(5,1);

    % Hepatic Cholesterol
    dxdt(1) = vBS * DietL ^ -0.1 * Cstor ^ -0.1 * Cpt ^ -0.1 + DietL +  v_pl2 * Cpla * Cpt ^ a4 - v_d * Chep *  DietL ^ 0.1  +  v_p * Cstor * Cpt ^ a3 * DietL ^ -0.5 - v_s * Chep - v_p14 * Chep;

    % Storage
    dxdt(2) = v_s * Chep - v_p * Cstor * Cpt ^ a3 * DietL ^ -0.5;

    % Peripheral Tissue Usage 
    dxdt(3) = v_pl1 * Cpla - v_pt * Cpt - v_ptt * Cpt;

    % Cholesterol Transport Plasma 
    dxdt(4) = v_pt * Cpt - v_pl1 * Cpla - v_pl2 * Cpla * Cpt ^ a4 + v_d * Chep * DietL ^ 0.5 - v_pl3 * Cpla * DietL ^ 0.5;

    % Estrogen Synthesis 
    dxdt(5) = v_pl3 * Cpla * DietL ^ 0.5 - v_es * Ces;
    
    % Return the derivatives
     dydt = dxdt;
end

