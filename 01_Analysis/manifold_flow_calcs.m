clear
clc
close all

% Inlet environment definition
Pc = 250 * 6894.76; % Chamber pressure, [Pa]
m_dot = 0.560640169 / 2; % Massflow fuel, [kg/s]
T_l = 293.16; % K

dP_pct_fu = 0.2; % Fuel stiffness
P_l = Pc * (1 + dP_pct_fu); % Fuel manifold pressure
rho_l = py.CoolProp.CoolProp.PropsSI('D','T', T_l,'P', P_l, 'ethanol');

% Manifold Geometry Definition
% A_max = 0.30653981 / 1550; % Ox [m^2]
% A_min = 0.21618775 / 1550;  % Ox [m^2]
A_max = 0.33381582 / 1550;  % Fu [m^2]
A_min = 0.23641132 / 1550; % Fu [m^2]

manifold_half_length = 1.65 * 0.0254 * pi; % [m]
e = 200e-6; % [micrometers]

% Manifold interpolation
steps = 100;
length = linspace(0, steps, steps);
manifold_mdot = m_dot - 0.9 * (m_dot - m_dot/steps) ./ steps .* (length); % mdot left after entrance through windows
cubic_area = (length + steps/2) .* (length - steps).^2 ./ (steps^2 * steps/2) * (A_max - A_min) + A_min;

% Velocity calculations
velo = manifold_mdot ./ (cubic_area * rho_l);

%% Pressure loss
deltax = manifold_half_length / steps; 
mu_lb = py.CoolProp.CoolProp.PropsSI('V','T', T_l, 'P', P_l, 'ethanol'); % viscosity of bulk coolant [Pa-s]
hydraulic_D = sqrt(cubic_area / pi);

Re_l = rho_l .* velo .* hydraulic_D ./ mu_lb;

P_l_manifold(1) = P_l;

for i = 2:1:steps
    colebrook = @(f) 1 / sqrt(f) + 2 * log10(e / (3.7 * hydraulic_D(i)) + 2.51 / (Re_l(i) * sqrt(f)));
    f = fzero(colebrook, 0.1);
    cf = f/4; % friction coefficient
    dP = 2 * cf * (deltax./hydraulic_D(i)) * rho_l * velo(i)^2;
    P_l_manifold(i) = P_l_manifold(i-1) - dP;
end
    

%% Figures

figure
hold on 
set(gcf,'color','w')
set(gca,'TickLabelInterpreter','latex')

plot(length, P_l_manifold(1) / 6895 - P_l_manifold / 6895, 'Linewidth', 2)

ylabel("dP [Psi]", 'Interpreter', 'latex')
xlabel("Length [m]", 'Interpreter', 'latex')
title("Manifold Pressure Loss",'Interpreter','latex')


figure
hold on 
set(gcf,'color','w')
set(gca,'TickLabelInterpreter','latex')

% plot(length, linspace(P_l / 6895, P_l / 6895, steps), 'Linewidth', 2)
plot(length, 1 / 2 * rho_l * velo .^ 2 / 6895, 'Linewidth', 2)

ylabel("Pressure [Psi]", 'Interpreter', 'latex')
xlabel("Length [m]", 'Interpreter', 'latex')
title("Dynamic Pressure",'Interpreter','latex')
% legend("Static", "Dynamic",'Interpreter','latex')



