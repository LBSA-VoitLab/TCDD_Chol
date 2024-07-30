% Unit: day
clear all
clc
%% ------------------------- DEFAULT PARAMETERS -----------------------------%%
%Hormone half-lives: https://www.kup.at/kup/pdf/4975.pdf
p.k1 = log(2)/(1/60);
p.k2 = log(2)/(1/60); %Assuming the conversion from Chol to P5 is fast such that the half-life of Chol is 1 minute.
p.k2a = 1;
p.k2b = 1;
p.k3 = log(2)/(1/60); %Assuming the conversion from P5 to P4 is fast such that the half-life of P5 is 1 minute.
p.k4 = 0.015 * log(2)/(1/60)*0.5; %Assuming the loss of P4 due to conversion to A4 and to secretion into circulation is fast such that the half-life of P4 is 1 minute, and the split is 50/50.
p.k5 = log(2)/(1/60)*0.5; 
p.k6 = log(2)/(1/60); %k6 and k7 are set this way so the difference between the two fluxes will be equal to k5 flux to maintain flux balance when A4 and T are both at 1.
p.k7 = log(2)/(1/60)*0.5; 
p.k8 = log(2)/(1/60)*0.5; 
p.k9 = log(2)/(1/60)*0.5;
p.k10 = 2*log(2)/2; % Estradiol half-life in circulation is 2 hours.
p.k11 = 0;
p.k12 = 0*log(2)/(1/60)*0.5; 
p.k13 = 0;
p.k14 = 0.125*log(2)/(5/60); %Progesterone half-life in circulation is about 5 min according to wikipedia

p.k15 = 2*1/1.45;
p.k17 = 2*log(2)/1; % LH half-life in circulation is 1 hour.

p.k18 = 0.1 * p.k15;
p.k19 = 10*log(2)/0.1;
p.k20 = log(2)/4;
p.k21 = 0.85*log(2)/5;
p.k22 = 2 * p.k9;
p.k23 = log(2)/1; % S half-life is 1 hour.
p.k24 = 1;
p.k25 = 1e4;
p.k26 = 3;
p.k27 = log(2)/0.1; % S1 half-life is 0.1 hour.
p.k28 = 40;
p.k29 = 2;
p.k31 = 0.25;
p.k30 = 10*p.k31;
p.k32 = 10*p.k20;
p.k33 = log(2)/(1/60);
p.k34 = 1/1.45;
p.k35 = 3.5*0.1 * p.k34;
p.k36 = log(2)/2;
p.k37 = 1e4;
p.k38 = 4;

p.k16 = 10 * p.k23;
p.k160 = p.k23;
p.k190 = 0.01*log(2)/0.1;


p.J16 = 15;
p.J19 = 16;
p.J22 = 3;
p.J25 = 5;
p.J29 = 5;
p.J30 = 1;
p.J32 = 9;
p.J34 = 3;
p.J37 = 5;
p.J38 = 5;


p.n16 = 5;
p.n19 = 5;
p.n22 = 100;
p.n25 = 100;
p.n29 = 10;
p.n30 = 5;
p.n32 = 10;
p.n34 = 1;
p.n37 = 100;
p.n38 = 1;

p.AFmax = 7.5;

p.Cholp = 1;

p.scale_E2p = 1/60; %Used to scale E2p to 1.
p.scale_P4p = 1/2.5; %Used to scale E2p to 1.

default_p = p;

%% ------------------------- INITIAL CONDITION-----------------------------%%
init.Chol= 1; %Cholesterol in ovary
init.P5  = 1; %Pregnenolone
init.P4  = 1; %Progesterone
init.A4  = 1; %Androstenedione
init.T   = 1; %Testosterone
init.E2  = 1; %Estradiol
init.E2p = 1.72;%2; %Estradiol in plasma
init.P4p = 0.62;%0; %Progesterone in plasma
init.LH  = 10; %LH in pituitary
init.LHp = 1; %LH in plasma
init.CL  = 0; %Corpus Luteum
init.AF  = 1; %Growing antral follicle
init.S   = 1; %Bistable signal S mediating E2p positive feedback onto LH release
init.S1   = 0; %Bistable signal S1 mediating LH surge-triggered AF collapse and CL formation
init.S2  = 0; %Intermediate promoting CL atresia
init.FSH  = 15; %FSH in pituitary
init.FSHp = 4; %FSH in plasma

default_init = init;

%% ------------------------- DEFAULT (AVERAGE MODEL) SIMULATION ---------------------------%%
p = default_p;
init = default_init;

y0 = cell2mat(struct2cell(init));
tspan = [0:0.1:28];
options = [];
[t,y] = ode15s('Steroidogenesis_model_ODE',tspan, y0, options, p);

Chol= y(:,1); %Cholesterol in ovary
P5  = y(:,2); %Pregnenolone
P4  = y(:,3); %Progesterone
A4  = y(:,4); %Androstenedione
T   = y(:,5); %Testosterone
E2  = y(:,6); %Estradiol
E2p = y(:,7); %Estradiol in plasma
P4p = y(:,8); %Progesterone in plasma
LH  = y(:,9); %LH in pituitary
LHp = y(:,10); %LH in plasma
CL  = y(:,11); %Corpus Luteum
AF  = y(:,12); %Growing antral follicle
S   = y(:,13); %Bistable signal S mediating E2p positive feedback onto LH release
S1  = y(:,14); %Bistable signal S1 mediating LH surge-triggered AF collapse and CL formation
S2  = y(:,15); %Intermediatepromoting CL atresia
FSH  = y(:,16); %FSH in pituitary
FSHp = y(:,17); %FSH in plasma

figure(100)
plot(t, E2p, 'r', 'LineWidth',2, 'DisplayName','E2') %Plasma estradiol
hold on
plot(t, P4p, 'b', 'LineWidth',2,  'DisplayName','P4') %Plasma progesterone
plot(t, LHp, 'g', 'LineWidth',2,  'DisplayName','LH') %Plasma LH
plot(t, FSHp, 'm', 'LineWidth',2,  'DisplayName','FSH') %Plasma FSH
%plot(t, AF, 'c', 'DisplayName','AF') %Antral follicle
%plot(t, CL, 'm', 'DisplayName','CL') %Corpus luteum
%plot(t, LH, 'm', 'DisplayName','LH') %Pituitary LH
%plot(t, FSH, 'm', 'DisplayName','FSH') %Pituitary FSH
%plot(t, S, 'k', 'DisplayName','S') %S
%plot(t, S1, 'k.', 'DisplayName','S1') %S
%plot(t, S2, 'k-', 'DisplayName','S2') %S
xlabel('Time (day)')
ylabel('Hormone level')
xticks([0:7:28])
legend

% figure(110)
% plot(t, LH, 'r') %Pituitary LH
% hold on
% xlabel('Time (day)')
% ylabel('Hormone level')
% legend('Pituitary LH')
% 
% figure(120)
% plot(t, Chol, 'r')
% hold on
% plot(t, P5, 'g')
% plot(t, P4, 'b')
% plot(t, A4, 'k')
% plot(t, T, 'c')
% plot(t, E2, 'm')
% xlabel('Time (day)')
% ylabel('Metabolite level')
% legend('Chol', 'P5', 'P4', 'A4', 'T', 'E2')
%     
%    figure(130)
%     plot(t, y(:,11), 'r') %corpus luteum
%     hold on
%     xlabel('Time (day)')
%     ylabel('corpus luteum')
%     legend('corpus luteum')

%Find the follicular phase length
[LHpmax,idx] = max(LHp)
LH_peak_time_default = t(idx)

%% ------------------------- Sensitivity Analysis for follicular phase length -----------%
percent_change = 0.1;
sensitivity_coeff_follicular_phase_length = []; 

%Run model at default parameter values
p = default_p;
init = default_init;

y0 = cell2mat(struct2cell(init));
tspan = [0:0.1:28];
options = [];
[t1,y1] = ode15s('Steroidogenesis_model_ODE',tspan, y0, options, p);
LHp1 = y1(:,10); %LH in plasma


%Find the follicular phase length
[LHpmax1,idx1] = max(LHp1)
LH_peak_time_default = t(idx1);


%-----------Run model at modified parameter values

%Compile parameter list for sensitivity analysis;
param_cell = struct2cell(default_p);
param_names = fieldnames(default_p);

for n = 1:length(param_cell)
    n
    
    %Increase the parameter value
        p = default_p;
        init = default_init;
        if n ~= (length(param_cell)+1)
            eval(strcat('p.', param_names{n}, '=', num2str(param_cell{n}), '*(1 + percent_change);'));
        end
        y0 = cell2mat(struct2cell(init));
        [t2,y2] = ode15s('Steroidogenesis_model_ODE',tspan, y0, options, p);
        LHp2 = y2(:,10); %LH in plasma
    
        %Find the follicular phase length
        [LHpmax2,idx2] = max(LHp2);
        LH_peak_time_up = t(idx2);
       
    
    %Decrease the parameter value
        p = default_p;
        init = default_init;
        if n ~= (length(param_cell)+1)
            eval(strcat('p.', param_names{n}, '=', num2str(param_cell{n}), '*(1 - percent_change);'));
        end
        y0 = cell2mat(struct2cell(init));
        [t3,y3] = ode15s('Steroidogenesis_model_ODE',tspan, y0, options, p);
        LHp3 = y3(:,10); %LH in plasma
    
        %Find the follicular phase length
        [LHpmax3,idx3] = max(LHp3)
        LH_peak_time_down = t(idx3); 
    
    %Calculate sensitivity coefficient
        sensitivity_coeff_follicular_phase_length(n) = mean([(LH_peak_time_up - LH_peak_time_default)/LH_peak_time_default/percent_change, (LH_peak_time_default - LH_peak_time_down)/LH_peak_time_default/percent_change])
end

%  tornado plot;
[~,idx] = sort(abs(sensitivity_coeff_follicular_phase_length), 'ascend'); % Obtain index after sorting by absolute values; 
param_names_follicular_phase_length = param_names(idx);
figure(1001)
barh(sensitivity_coeff_follicular_phase_length(idx));
ylim([0.5,length(param_cell)+1.5]);
xlabel('sensitivity coefficient (follicular phase length)');
yticks([1:1:length(param_cell)+1])
yticklabels(param_names_follicular_phase_length);
pbaspect([0.75 1 1])

%% ------------------------- MONTE CARLO POPULATION SIMULATION ---------------------------%%
num_indv    = 100; %100 women

p = default_p;
init = default_init;

init_E2p_out = [];
init_P4p_out = [];
init_LHp_out = [];
init_FSHp_out = [];
J16_out = [];
k10_out = [];
k14_out = [];
k21_out = [];
k18_out = [];

out   = {};
for n = 1:1:num_indv
    n
    
    
    %Initial E2p
    m = default_init.E2p; %Mean
    v = (m*0.1)^2; %variance
    E2p_mean = log((m^2)/sqrt(v+m^2));
    E2p_sd = sqrt(log(v/(m^2)+1)); 
    init.E2p = lognrnd(E2p_mean, E2p_sd, 1);
    init_E2p_out = [init_E2p_out; init.E2p]; 
%     
%     %Initial P4p
%     m = default_init.P4p; %Mean
%     v = (m*0.1)^2; %variance
%     P4p_mean = log((m^2)/sqrt(v+m^2));
%     P4p_sd = sqrt(log(v/(m^2)+1)); 
%     init.P4p = lognrnd(P4p_mean, P4p_sd, 1);
%     init_P4p_out = [init_P4p_out; init.P4p]; 
    
    %Initial LHp
    m = default_init.LHp; %Mean
    v = (m*0.1)^2; %variance
    LHp_mean = log((m^2)/sqrt(v+m^2));
    LHp_sd = sqrt(log(v/(m^2)+1)); 
    init.LHp = lognrnd(LHp_mean, LHp_sd, 1);
    init_LHp_out = [init_LHp_out; init.LHp]; 
    
    %Initial FSHp
    m = default_init.FSHp; %Mean
    v = (m*0.1)^2; %variance
    FSHp_mean = log((m^2)/sqrt(v+m^2));
    FSHp_sd = sqrt(log(v/(m^2)+1)); 
    init.FSHp = lognrnd(FSHp_mean, FSHp_sd, 1);
    init_FSHp_out = [init_FSHp_out; init.FSHp]; 
        
    
    %k10 sampling for plasma E2 clearance
    m	= default_p.k10; %Mean
    v = (m*0.05)^2; %variance 0.1
    k10_mean = log((m^2)/sqrt(v+m^2));
    k10_sd = sqrt(log(v/(m^2)+1));
    p.k10 = lognrnd(k10_mean, k10_sd, 1);
    k10_out = [k10_out; p.k10];
    
    %k14 sampling for plasma P4 clearance
    m	= default_p.k14; %Mean
    v = (m*0.05)^2; %variance 0.1
    k14_mean = log((m^2)/sqrt(v+m^2));
    k14_sd = sqrt(log(v/(m^2)+1));
    p.k14 = lognrnd(k14_mean, k14_sd, 1);
    k14_out = [k14_out; p.k14];
        
    %k21 sampling for AF growth rate
    m	= default_p.k21; %Mean
    v = (m*0.05)^2; %variance 0.1
    k21_mean = log((m^2)/sqrt(v+m^2));
    k21_sd = sqrt(log(v/(m^2)+1));
    p.k21 = lognrnd(k21_mean, k21_sd, 1);
    k21_out = [k21_out; p.k21];
    
    
    %k18 sampling for different amount of pituitary LH stored such that the LH surge can reach different peak levels
    m	= default_p.k18; %Mean
    v = (m*0.2)^2; %variance
    k18_mean = log((m^2)/sqrt(v+m^2));
    k18_sd = sqrt(log(v/(m^2)+1));
    p.k18 = lognrnd(k18_mean, k18_sd, 1);
    k18_out = [k18_out; p.k18];
    
    %J16 sampling for E2 threshold triggering LH surge
    m	= default_p.J16; %Mean
    v = (m*0.05)^2; %variance
    J16_mean = log((m^2)/sqrt(v+m^2));
    J16_sd = sqrt(log(v/(m^2)+1));
    p.J16 = lognrnd(J16_mean, J16_sd, 1);
    J16_out = [J16_out; p.J16];
    
    
    
    y0 = cell2mat(struct2cell(init));
    tspan = [0:0.1:28];
    options = [];
    [t,y] = ode15s('Steroidogenesis_model_ODE',tspan, y0, options, p);
    
    Chol= y(:,1); %Cholesterol in ovary
    P5  = y(:,2); %Pregnenolone
    P4  = y(:,3); %Progesterone
    A4  = y(:,4); %Androstenedione
    T   = y(:,5); %Testosterone
    E2  = y(:,6); %Estradiol
    E2p = y(:,7); %Estradiol in plasma
    P4p = y(:,8); %Progesterone in plasma
    LH  = y(:,9); %LH in pituitary
    LHp = y(:,10); %LH in plasma
    CL  = y(:,11); %Corpus Luteum
    AF  = y(:,12); %Growing antral follicle
    S   = y(:,13); %Bistable signal S mediating E2p positive feedback onto LH release
    S1  = y(:,14); %Bistable signal S1 mediating LH surge-triggered AF collapse and CL formation
    S2  = y(:,15); %Intermediatepromoting CL atresia
    FSH  = y(:,16); %FSH in pituitary
    FSHp = y(:,17); %FSH in plasma
    
    out.P4p(:,n) = P4p;
    out.E2p(:,n) = E2p;
    out.LHp(:,n) = LHp;
    out.FSHp(:,n) = FSHp;
    
%     %Plasma E2
%     figure(199)
%     plot(t, E2p, 'DisplayName','Plasma E2') %Plasma estradiol
%     hold on
%     %plot(t, P4p, 'b', 'DisplayName','P4p') %Plasma progesterone
%     %plot(t, LHp, 'g', 'DisplayName','LHp') 
%     %plot(t, AF, 'c', 'DisplayName','AF') 
%     %plot(t, CL, 'm', 'DisplayName','CL') 
%     plot(t, LH, 'm', 'DisplayName','LH') 
%     %plot(t, S, 'k', 'DisplayName','S') 
%     %plot(t, S1, 'k.', 'DisplayName','S1') 
%     %plot(t, S2, 'k-', 'DisplayName','S2') 
%     xlabel('Time (day)')
%     ylabel('Hormone level')
%     legend
%     
%     %Plasma P4
%     figure(210)
%     plot(t, P4p, 'DisplayName','Plasma P4') 
%     hold on
%     xlabel('Time (day)')
%     ylabel('Hormone level')
%     legend
%     
%     %Plasma LH
%     figure(220)
%     plot(t, LHp, 'DisplayName','Plasma LH')
%     hold on
%     xlabel('Time (day)')
%     ylabel('Hormone level')
%     legend
    
end
    
%Time-course figures
%Plasam E2
figure(200)
plot(t,out.E2p(:,1:end), 'Color', [0.7 0.7 0.7])
hold on
xlabel('Day')
ylabel('Plasam E2')
percentile_E2p = prctile(out.E2p',[2.5,50,97.5])';
%plot(t,percentile_E2p(:,2), 'Color', 'blue' ,'LineWidth',3)
plot(t,percentile_E2p(:,[1,3]), 'b--', 'LineWidth',3)
plot(t,mean(out.E2p'), 'Color', 'blue' ,'LineWidth',3)
title('Plasma E2')
xlim([0 28])

%Plasam P4
figure(210)
plot(t,out.P4p(:,1:end), 'Color', [0.7 0.7 0.7])
hold on
xlabel('Day')
ylabel('Plasam P4')
percentile_P4p = prctile(out.P4p',[2.5,50,97.5])';
%plot(t,percentile_P4p(:,2), 'Color', 'blue' ,'LineWidth',3)
plot(t,percentile_P4p(:,[1,3]), 'b--', 'LineWidth',3)
plot(t,mean(out.P4p'), 'Color', 'blue' ,'LineWidth',3)
title('Plasma P4')
xlim([0 28])

%Plasam LH
figure(220)
plot(t,out.LHp(:,1:end), 'Color', [0.7 0.7 0.7])
hold on
xlabel('Day)')
ylabel('Plasma LH')
percentile_LHp = prctile(out.LHp',[2.5,50,97.5])';
%plot(t,percentile_LHp(:,2), 'Color', 'blue' ,'LineWidth',3)
plot(t,percentile_LHp(:,[1,3]), 'b--', 'LineWidth',3)
plot(t,mean(out.LHp'), 'Color', 'blue' ,'LineWidth',3)
title('Plasma LH')
xlim([0 28])

%Plasam FSH
figure(230)
plot(t,out.FSHp(:,1:end), 'Color', [0.7 0.7 0.7])
hold on
xlabel('Day)')
ylabel('Plasma LH')
percentile_FSHp = prctile(out.FSHp',[2.5,50,97.5])';
%plot(t,percentile_FSHp(:,2), 'Color', 'blue' ,'LineWidth',3)
plot(t,percentile_FSHp(:,[1,3]), 'b--', 'LineWidth',3)
plot(t,mean(out.FSHp'), 'Color', 'blue' ,'LineWidth',3)
title('Plasma FSH')
xlim([0 28])

%Capture the follicular phase length
LH_peak_time =[];
for n = 1:1:num_indv
    %Find the LH time
    [LHpmax,idx] = max(out.LHp(:,n))
    LH_peak_time = [LH_peak_time, t(idx)]; %code is breaking so run outside of loop or change tspan
end
figure(240) 
h = histogram(LH_peak_time,'Normalization','probability','BinWidth',0.5);
hold on
%nbins = h.NumBins to find number of bins
probability = h.Values;
days = h.BinEdges(1:end-1);  
xlabel('Folliclular phase length (day)')
ylabel('Probability')
mean_folliclular_phase_length   = mean(LH_peak_time);
sd_folliclular_phase_length    = std(LH_peak_time);
percentile_folliclular_phase_length = prctile(LH_peak_time,[2.5,50,97.5]);
median_folliclular_phase_length   = median(LH_peak_time);
title('Distribution of Folliclular Phase Length')
text(16, 0.14, strcat('mean=', num2str(round(mean_folliclular_phase_length,1)), '+/-', num2str(round(sd_folliclular_phase_length,1))));

%% ------------------------- DIOXIN EFFECT ON AVERAGE MODEL ---------------------------%%
p = default_p;
init = default_init;

p.k2 = default_p.k2 * 0.8;

%p.k6 = default_p.k6 * 0.01;

y0 = cell2mat(struct2cell(init));
tspan = [0:0.1:28];
options = [];
[t,y] = ode15s('Steroidogenesis_model_ODE',tspan, y0, options, p);

Chol= y(:,1); %Cholesterol in ovary
P5  = y(:,2); %Pregnenolone
P4  = y(:,3); %Progesterone
A4  = y(:,4); %Androstenedione
T   = y(:,5); %Testosterone
E2  = y(:,6); %Estradiol
E2p = y(:,7); %Estradiol in plasma
P4p = y(:,8); %Progesterone in plasma
LH  = y(:,9); %LH in pituitary
LHp = y(:,10); %LH in plasma
CL  = y(:,11); %Corpus Luteum
AF  = y(:,12); %Growing antral follicle
S   = y(:,13); %Bistable signal S mediating E2p positive feedback onto LH release
S1  = y(:,14); %Bistable signal S1 mediating LH surge-triggered AF collapse and CL formation
S2  = y(:,15); %Intermediatepromoting CL atresia
FSH  = y(:,16); %FSH in pituitary
FSHp = y(:,17); %FSH in plasma

figure(301)
plot(t, E2p, 'r', 'LineWidth',2, 'DisplayName','E2') %Plasma estradiol
hold on
plot(t, P4p, 'b', 'LineWidth',2,  'DisplayName','P4') %Plasma progesterone
plot(t, LHp, 'g', 'LineWidth',2,  'DisplayName','LH') %Plasma LH
plot(t, FSHp, 'm', 'LineWidth',2,  'DisplayName','FSH') %Plasma FSH
%plot(t, AF, 'c', 'DisplayName','AF') %Antral follicle
%plot(t, CL, 'm', 'DisplayName','CL') %Corpus luteum
%plot(t, LH, 'm', 'DisplayName','LH') %Pituitary LH
%plot(t, FSH, 'm', 'DisplayName','FSH') %Pituitary FSH
%plot(t, S, 'k', 'DisplayName','S') %S
%plot(t, S1, 'k.', 'DisplayName','S1') %S
%plot(t, S2, 'k-', 'DisplayName','S2') %S
xlabel('Time (day)')
ylabel('Hormone level')
xticks([0:7:28])
legend

% figure(110)
% plot(t, LH, 'r') %Pituitary LH
% hold on
% xlabel('Time (day)')
% ylabel('Hormone level')
% legend('Pituitary LH')
% 
% figure(120)
% plot(t, Chol, 'r')
% hold on
% plot(t, P5, 'g')
% plot(t, P4, 'b')
% plot(t, A4, 'k')
% plot(t, T, 'c')
% plot(t, E2, 'm')
% xlabel('Time (day)')
% ylabel('Metabolite level')
% legend('Chol', 'P5', 'P4', 'A4', 'T', 'E2')
%     
%    figure(130)
%     plot(t, y(:,11), 'r') %corpus luteum
%     hold on
%     xlabel('Time (day)')
%     ylabel('corpus luteum')
%     legend('corpus luteum')

%Find the follicular phase length
[LHpmax,idx] = max(LHp)
LH_peak_time_dioxin = t(idx)

%% ------------------------- DIOXIN EFFECT ON MONTE CARLO POPULATION SIMULATION ---------------------------%%
num_indv    = 200; %100 women

p = default_p;
init = default_init;

TCDD_vector = [1,10,100,1000,10000];
a= 0.77
Ki = 100 %ppt
n = 0.5
mean_folliclular_phase_length_vector = [];
sd_folliclular_phase_length_vector = [];
for i = 1:length(TCDD_vector)
    i
    TCDD = TCDD_vector(i);
    k2_fold_change = a + (1-a) * Ki^n / (Ki^n + TCDD^n)
    p.k2 = default_p.k2 * k2_fold_change;


    init_E2p_out = [];
    init_P4p_out = [];
    init_LHp_out = [];
    init_FSHp_out = [];
    J16_out = [];
    k10_out = [];
    k14_out = [];
    k21_out = [];
    k18_out = [];

    out   = {};
    for j = 1:1:num_indv
        j


        %Initial E2p
        m = default_init.E2p; %Mean
        v = (m*0.1)^2; %variance
        E2p_mean = log((m^2)/sqrt(v+m^2));
        E2p_sd = sqrt(log(v/(m^2)+1)); 
        init.E2p = lognrnd(E2p_mean, E2p_sd, 1);
        init_E2p_out = [init_E2p_out; init.E2p]; 
    %     
    %     %Initial P4p
    %     m = default_init.P4p; %Mean
    %     v = (m*0.1)^2; %variance
    %     P4p_mean = log((m^2)/sqrt(v+m^2));
    %     P4p_sd = sqrt(log(v/(m^2)+1)); 
    %     init.P4p = lognrnd(P4p_mean, P4p_sd, 1);
    %     init_P4p_out = [init_P4p_out; init.P4p]; 

        %Initial LHp
        m = default_init.LHp; %Mean
        v = (m*0.1)^2; %variance
        LHp_mean = log((m^2)/sqrt(v+m^2));
        LHp_sd = sqrt(log(v/(m^2)+1)); 
        init.LHp = lognrnd(LHp_mean, LHp_sd, 1);
        init_LHp_out = [init_LHp_out; init.LHp]; 

        %Initial FSHp
        m = default_init.FSHp; %Mean
        v = (m*0.1)^2; %variance
        FSHp_mean = log((m^2)/sqrt(v+m^2));
        FSHp_sd = sqrt(log(v/(m^2)+1)); 
        init.FSHp = lognrnd(FSHp_mean, FSHp_sd, 1);
        init_FSHp_out = [init_FSHp_out; init.FSHp]; 


        %k10 sampling for plasma E2 clearance
        m	= default_p.k10; %Mean
        v = (m*0.05)^2; %variance 0.1
        k10_mean = log((m^2)/sqrt(v+m^2));
        k10_sd = sqrt(log(v/(m^2)+1));
        p.k10 = lognrnd(k10_mean, k10_sd, 1);
        k10_out = [k10_out; p.k10];

        %k14 sampling for plasma P4 clearance
        m	= default_p.k14; %Mean
        v = (m*0.05)^2; %variance 0.1
        k14_mean = log((m^2)/sqrt(v+m^2));
        k14_sd = sqrt(log(v/(m^2)+1));
        p.k14 = lognrnd(k14_mean, k14_sd, 1);
        k14_out = [k14_out; p.k14];

        %k21 sampling for AF growth rate
        m	= default_p.k21; %Mean
        v = (m*0.05)^2; %variance 0.1
        k21_mean = log((m^2)/sqrt(v+m^2));
        k21_sd = sqrt(log(v/(m^2)+1));
        p.k21 = lognrnd(k21_mean, k21_sd, 1);
        k21_out = [k21_out; p.k21];


        %k18 sampling for different amount of pituitary LH stored such that the LH surge can reach different peak levels
        m	= default_p.k18; %Mean
        v = (m*0.2)^2; %variance
        k18_mean = log((m^2)/sqrt(v+m^2));
        k18_sd = sqrt(log(v/(m^2)+1));
        p.k18 = lognrnd(k18_mean, k18_sd, 1);
        k18_out = [k18_out; p.k18];

        %J16 sampling for E2 threshold triggering LH surge
        m	= default_p.J16; %Mean
        v = (m*0.05)^2; %variance
        J16_mean = log((m^2)/sqrt(v+m^2));
        J16_sd = sqrt(log(v/(m^2)+1));
        p.J16 = lognrnd(J16_mean, J16_sd, 1);
        J16_out = [J16_out; p.J16];



        y0 = cell2mat(struct2cell(init));
        tspan = [0:0.1:28];
        options = [];
        [t,y] = ode15s('Steroidogenesis_model_ODE',tspan, y0, options, p);

        Chol= y(:,1); %Cholesterol in ovary
        P5  = y(:,2); %Pregnenolone
        P4  = y(:,3); %Progesterone
        A4  = y(:,4); %Androstenedione
        T   = y(:,5); %Testosterone
        E2  = y(:,6); %Estradiol
        E2p = y(:,7); %Estradiol in plasma
        P4p = y(:,8); %Progesterone in plasma
        LH  = y(:,9); %LH in pituitary
        LHp = y(:,10); %LH in plasma
        CL  = y(:,11); %Corpus Luteum
        AF  = y(:,12); %Growing antral follicle
        S   = y(:,13); %Bistable signal S mediating E2p positive feedback onto LH release
        S1  = y(:,14); %Bistable signal S1 mediating LH surge-triggered AF collapse and CL formation
        S2  = y(:,15); %Intermediatepromoting CL atresia
        FSH  = y(:,16); %FSH in pituitary
        FSHp = y(:,17); %FSH in plasma

        out.P4p(:,j) = P4p;
        out.E2p(:,j) = E2p;
        out.LHp(:,j) = LHp;

    %     %Plasma E2
    %     figure(200)
    %     plot(t, E2p, 'DisplayName','Plasma E2') %Plasma estradiol
    %     hold on
    %     %plot(t, P4p, 'b', 'DisplayName','P4p') %Plasma progesterone
    %     %plot(t, LHp, 'g', 'DisplayName','LHp') %Plasma LH
    %     %plot(t, AF, 'c', 'DisplayName','AF') %Plasma LH
    %     %plot(t, CL, 'm', 'DisplayName','CL') %Plasma LH
    %     %plot(t, LH, 'm', 'DisplayName','LH') %Plasma LH
    %     %plot(t, S, 'k', 'DisplayName','S') %S
    %     %plot(t, S1, 'k.', 'DisplayName','S1') %S
    %     %plot(t, S2, 'k-', 'DisplayName','S2') %S
    %     xlabel('Time (day)')
    %     ylabel('Hormone level')
    %     legend
    %     
    %     %Plasma P4
    %     figure(210)
    %     plot(t, P4p, 'DisplayName','Plasma P4') 
    %     hold on
    %     xlabel('Time (day)')
    %     ylabel('Hormone level')
    %     legend
    %     
    %     %Plasma LH
    %     figure(220)
    %     plot(t, LHp, 'DisplayName','Plasma LH')
    %     hold on
    %     xlabel('Time (day)')
    %     ylabel('Hormone level')
    %     legend

    end

    %Time-course figures
    %Plasam E2
    figure(400)
    plot(t,out.E2p(:,1:end), 'Color', [0.7 0.7 0.7])
    hold on
    xlabel('Day')
    ylabel('Plasam E2')
    percentile_E2p = prctile(out.E2p',[2.5,50,97.5])';
    %plot(t,percentile_E2p(:,2), 'Color', 'blue' ,'LineWidth',3)
    plot(t,percentile_E2p(:,[1,3]), 'b--', 'LineWidth',3)
    plot(t,mean(out.E2p'), 'Color', 'blue' ,'LineWidth',3)
    title('Plasma E2')
    xlim([0 28])

    %Plasam P4
    figure(410)
    plot(t,out.P4p(:,1:end), 'Color', [0.7 0.7 0.7])
    hold on
    xlabel('Day')
    ylabel('Plasam P4')
    percentile_P4p = prctile(out.P4p',[2.5,50,97.5])';
    %plot(t,percentile_P4p(:,2), 'Color', 'blue' ,'LineWidth',3)
    plot(t,percentile_P4p(:,[1,3]), 'b--', 'LineWidth',3)
    plot(t,mean(out.P4p'), 'Color', 'blue' ,'LineWidth',3)
    title('Plasma P4')
    xlim([0 28])

    %Plasam LH
    figure(420)
    plot(t,out.LHp(:,1:end), 'Color', [0.7 0.7 0.7])
    hold on
    xlabel('Day)')
    ylabel('Plasma LH')
    percentile_LHp = prctile(out.LHp',[2.5,50,97.5])';
    %plot(t,percentile_LHp(:,2), 'Color', 'blue' ,'LineWidth',3)
    plot(t,percentile_LHp(:,[1,3]), 'b--', 'LineWidth',3)
    plot(t,mean(out.LHp'), 'Color', 'blue' ,'LineWidth',3)
    title('Plasma LH')
    xlim([0 28])

    %Capture the follicular phase length
    LH_peak_time =[];
    for j = 1:1:num_indv
        %Find the LH time
        [LHpmax,idx] = max(out.LHp(:,j))
        LH_peak_time = [LH_peak_time, t(idx)]; %code is breaking so run outside of loop or change tspan
    end
    figure(430) 
    h = histogram(LH_peak_time,'Normalization','probability','BinWidth',0.5);
    hold on
    %nbins = h.NumBins to find number of bins
    probability = h.Values;
    days = h.BinEdges(1:end-1);
    xlabel('Folliclular phase length (day)')
    ylabel('Probability')
    mean_folliclular_phase_length   = mean(LH_peak_time);
    mean_folliclular_phase_length_vector = [mean_folliclular_phase_length_vector, mean_folliclular_phase_length];
    sd_folliclular_phase_length    = std(LH_peak_time);
    sd_folliclular_phase_length_vector = [sd_folliclular_phase_length_vector, sd_folliclular_phase_length];
    percentile_folliclular_phase_length = prctile(LH_peak_time,[2.5,50,97.5]);
    median_folliclular_phase_length   = median(LH_peak_time);
    title('Cumulative Distribution of Folliclular Phase Length')
    text(16, 0.14, strcat('mean=', num2str(round(mean_folliclular_phase_length,1)), '+/-', num2str(round(sd_folliclular_phase_length,1))));
    pbaspect([1.5, 1, 1])
    
    %Cumulative probability plot
    figure(440)
    hold on
    h = cdfplot(LH_peak_time)
    pbaspect([1.5, 1, 1])

end

figure(450)
errorbar(TCDD_vector, mean_folliclular_phase_length_vector, sd_folliclular_phase_length_vector)
xlabel('TCDD serum level (ppt)')
ylabel('Follicular Phase Length (mean +/- sd)')