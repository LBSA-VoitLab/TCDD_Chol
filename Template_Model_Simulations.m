% Main script to call the ODE solver
clear

% Parameters
Par = [20955300; 1506; 6667 * 0.2; 77.32 * 0.2; 18460; 53200000; 259.81; 0.0000248 * 0.2; 2.5425 * 1.8; 10.33 * 1.8; 0.72 * 0.2; 517.26 / 1.93; -0.1; 0.5];

% Initial conditions
x0 = [18150; 0.56; 266.026; 3266.373; 517.26 / 1.93];

% Time span
tspan = 1:10000; 

% Solve the ODEs
tic
tspan = 1:10000; 
[t,x] = ode23s(@Template_model_ODES,tspan,x0',[],Par); 
toc %26s



% Plot for Cholesterol Biosynthesis
figure(1)
plot(t, x(:, 1))
xlabel('time')
ylabel('Cholesterol Biosynthesis')

% Plot for Storage
figure(2)
plot(t, x(:, 2))
xlabel('time')
ylabel('Storage')

% Plot for Peripheral Tissue Usage 
figure(3)
plot(t, x(:, 3))
xlabel('time')
ylabel('Peripheral Tissue Usage')

% Plot for Cholesterol Transport Plasma
figure(4)
plot(t, x(:, 4))
xlabel('time')
ylabel('Cholesterol Transport Plasma')

% Plot for Estrogen Synthesis 
figure(5)
plot(t, x(:, 5))
xlabel('time')
ylabel('Estrogen Synthesis')



% Report the steady-state values
steady_state_values = x(end, :);
disp('Steady-state values:')
disp('Hepatic Cholesterol:')
disp(steady_state_values(1))
disp('Storage:')
disp(steady_state_values(2))
disp('Peripheral Tissue Usage:')
disp(steady_state_values(3))
disp('Cholesterol Transport Plasma:')
disp(steady_state_values(4))
disp('Estrogen Synthesis:')
disp(steady_state_values(5))
