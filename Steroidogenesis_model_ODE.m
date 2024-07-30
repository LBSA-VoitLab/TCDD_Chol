function dydt = Steroidogenesis_model_ODE(t, y, option, p)


    Chol= y(1); %Cholesterol in ovary
    P5  = y(2); %Pregnenolone
    P4  = y(3); %Progesterone
    A4  = y(4); %Androstenedione
    T   = y(5); %Testosterone
    E2  = y(6); %Estradiol
    E2p = y(7); %Estradiol in plasma
    P4p = y(8); %Progesterone in plasma
    LH  = y(9); %LH in pituitary
    LHp = y(10); %LH in plasma
    CL  = y(11); %Corpus Luteum
    AF  = y(12); %Growing antral follicles
    S   = y(13); %Bistable signal S mediating E2p positive feedback onto LH release
    S1  = y(14); %Bistable signal S1 mediating LH surge-triggered AF collapse and CL formation
    S2  = y(15); %Intermediate promoting CL atresia
    FSH  = y(16); %FSH in pituitary
    FSHp = y(17); %FSH in plasma
    


  % Define the system of ODEs
    dydt = zeros(length(y),1);
    
    %Chol
    %dydt(1) = p.k1 * AF * p.Cholp + p.k33 * CL * p.Cholp - p.k2 * Chol;
    dydt(1) = 0; %p.k1 * p.Cholp + p.k33 * p.Cholp - p.k2 * Chol;
    
    %P5
    %dydt(2) = p.k2 * Chol - p.k3 * P5 - p.k11 * P5;
    dydt(2) = p.k2 * (p.k2a*AF + p.k2b*CL) * Chol - p.k3 * P5 - p.k11 * P5;
    
    %P4
    dydt(3) = p.k3 * P5 - p.k4 * CL * P4 - p.k5 * P4;
    
    %A4
    dydt(4) = p.k5 * P4 - p.k6 * A4 + p.k7 * T - p.k12 * A4;
    
    %T
    dydt(5) = p.k6 * A4 - p.k7 * T - p.k8 * T - p.k13 * T;
    
    %E2
    dydt(6) = p.k8 * T - p.k9 * E2 - p.k22 * LHp^p.n22/(p.J22^p.n22+LHp^p.n22) * E2;
    
    %E2p
    dydt(7) = p.scale_E2p * p.k9 * E2 - p.k10 * E2p;
    
    %P4p
    dydt(8) = p.scale_P4p * p.k4 * CL * P4 - p.k14 * P4p;
    
    %LH
    dydt(9) = p.k15 - p.k25 * S^p.n25/(p.J25^p.n25+S^p.n25) * LH - p.k18 * LH;
    
    %LHp
    dydt(10) = p.k25 * S^p.n25/(p.J25^p.n25+S^p.n25) * LH + p.k18 * LH - p.k17 * LHp;
    
    %CL
    dydt(11) = p.k29 * S1^p.n29/(p.J29^p.n29+S1^p.n29) - p.k20 * CL - p.k32 * S2^p.n32 / (p.J32^p.n32 + S2^p.n32) * CL;
    
    %AF
    dydt(12) = p.k21 * AF^2 * (1 - AF/p.AFmax) - p.k29 * S1^p.n29/(p.J29^p.n29+S1^p.n29) * AF; %Logistic function to describe exponential growth when AF is small and stop growing when AF is large 

    %S
    dydt(13) = p.k160 + p.k16 * S^p.n16 / ((p.J16/(p.k24+E2p))^p.n16 +  S^p.n16 ) - p.k23*S;
    
    %S1
    dydt(14) = p.k190 + p.k19 * S1^p.n19 / ((p.J19/(p.k26+p.k28*LHp))^p.n19 +  S1^p.n19 ) - p.k27*S1;
    
    %S2
    dydt(15) = p.k30 * S1^p.n30 / (p.J30^p.n30 + S1^p.n30) - p.k31*S2;
    
    %FSH
    %dydt(16) = p.k34 - p.k37 * S^p.n37/(p.J37^p.n37+S^p.n37) * FSH - p.k35 * FSH;
    dydt(16) = p.k34 * p.J34^p.n34/(p.J34^p.n34+P4p^p.n38) * p.k38 * p.J38^p.n38/(p.J38^p.n38+E2p^p.n38) - p.k37 * S^p.n37/(p.J37^p.n37+S^p.n37) * FSH - p.k35 * FSH;
    
    %FSHp
    dydt(17) = p.k37 * S^p.n37/(p.J37^p.n37+S^p.n37) * FSH + p.k35 * FSH - p.k36 * FSHp;
    
    
    
end