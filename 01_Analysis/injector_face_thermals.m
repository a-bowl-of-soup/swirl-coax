clear;
clc;
close all

Pc = 250 * 6894.76; % Chamber pressure, [Pa]
Pe = 17; % exit pressure [psi]
fuel_cea = 'C3H8O,2propanol';
T_l_cea = 293.15;
fuel = 'ethanol';
T_l = 293.15;
oxidizer = 'O2(L)';
oxidizer_temp = 90.17;
OF = 1.2;

% c_star, cf, isp, exp_ratio, M, gamma, P, T, rho, mu, Pr, Mw, k, son, cp
[~, ~, ~, ~, ~, gamma, ~, Tc_ns, ~, ~, Pr_g, ~, ~, ~, ~] = RunCEA(Pc / 6894.76, Pe, fuel_cea, 0, T_l_cea, oxidizer, oxidizer_temp, OF, 0, 0, 1, 0, 0, 'cea');

% Mass flow
manifold_mdot = 0.560640169; % Massflow fuel, [kg/s] 

% Injector definition  
t_w = 0.1 * 0.0254;  % [m]
e = 200e-6; % [micrometers]
dP_pct_fu = 0.2; % Fuel stiffness
P_l = Pc * (1 + dP_pct_fu)*2; % Fuel manifold pressure

% wall material properties
k_w = 25; % Thermal Conductivity of Wall [W/m-K]

% Reference diameters
rho_fu = py.CoolProp.CoolProp.PropsSI('D','T', T_l,'P', P_l, fuel); % Fu density in manifold, [kg/m^3]
P_dyn = P_l * 0.01;
velo = 1.3; %sqrt(2 * P_dyn / rho_fu); % [m/s]
hydraulic_D = sqrt(manifold_mdot / (velo * rho_fu) / pi) * 2;

% Iteration initialization
T_wg = 1000; % initial guess of wall side temperature [K]
steps = 1;
htc_array = linspace(2000, 2000, steps); % htc, [W/m^2-K]
M = 0.1;

qdot_tolerance = 0.0001;
for i = 1:steps
    converged = 0;
    T_wg_mn = 500; % minimum temperature bound
    T_wg_mx = 3000; % maximum temperature bound

    while ~(converged)
        h_g = htc_array(i);
        r = Pr_g ^ (1 / 3); % recovery factor 
        T_r = Tc_ns * (1 + (gamma - 1) / 2 * r * M ^ 2); % recovery temperature [K]
        qdot_g = h_g * (T_r - T_wg); % gas convective heat flux [W/m^2] 
        
        % Calculate liquid wall temperature
        T_wl = T_wg - qdot_g * t_w / k_w; % liquid wall temperature [K] 
        T_avg = (T_l + T_wl) / 2;

        % Calculate liquid film coefficient
        mu_lb = py.CoolProp.CoolProp.PropsSI('V','T', T_avg, 'P', P_l, fuel); % viscosity of bulk coolant [Pa-s]
        cp_l = py.CoolProp.CoolProp.PropsSI('C' , 'T', T_avg, 'P', P_l, fuel); % specific heat of coolant [J/kg-k] 
        k_l = py.CoolProp.CoolProp.PropsSI('L', 'T', T_avg, 'P', P_l, fuel); % thermal conductivity of coolant [W/m-K]
        
        Re_l = (4 * manifold_mdot) / (pi * hydraulic_D * mu_lb); % reynolds number for channel flow
        Pr_l = (cp_l * mu_lb) / k_l; % prantl number
        
        colebrook = @(f) 1 / sqrt(f) + 2 * log10(e / (3.7 * hydraulic_D) + 2.51 / (Re_l * sqrt(f)));
        f = fzero(colebrook, 0.1);
        Nu_l = (f / 8) * (Re_l - 1000) * Pr_l / (1 + 12.7 * (f / 8) ^ 0.5 * (Pr_l ^ (2 / 3) - 1)); % Gnielinksy correlation nusselt number [N/A] - 0.5 < Pr < 2000, 3000 < Re < 5e6        h_l = (Nu_l * k_l) / hydraulic_D; % liquid film coefficient [W/m^2-K] ALEX CITE SOURCE
        h_l(i) = (Nu_l * k_l) / hydraulic_D;
    
        % Step 8: Calculate liquid-side convective heat flux
        qdot_l = h_l(i) * (T_wl - T_l); % liquid convective heat flux [W/m^2] (Heister EQ 6.29).
        
        % Step 9: Check for convergence and continue loop / next step
        if abs(qdot_g - qdot_l) > qdot_tolerance % check for tolerance
            % convergence loop
            if qdot_g - qdot_l > 0
                T_wg_mn = T_wg;
            else 
                T_wg_mx = T_wg;
            end 
            T_wg = (T_wg_mx + T_wg_mn) / 2;
        else
            converged = 1;
            T_faceplate(i) = T_wg; 
            T_faceplate_cold(i) = T_wl; 
        end
    end
end
T_faceplate
%% Stress
E = 30e+9; 
CTE = 15e-6; 
v = 0.3; 
beta = 1; 
YTS = 105e6; 
e = 0.5;
ft_ult = 200e6; 

P = P_l - Pc;
r = hydraulic_D / 2;
dT1 = T_faceplate - T_faceplate_cold;
dT2 = (T_faceplate + T_faceplate_cold) / 2 - T_l; 

s_hoop = P_l * r / t_w; % Hoop stress, [Pa]
s_th = E * CTE * dT1 / (2 * (1 - v)); % Thermal stress in an unconstrained thin-walled cylinder, [Pa]

s_a = s_th + E .* CTE .* dT2 + P_l * r / (2 * t_w); % Stress axial (includes axial stress from pressure Pr/2t, you can delete that if it's not applicable), [Pa]
s_t = s_th * beta + s_hoop; % Stress tangential, [Pa]

s_eqv = sqrt(s_t .^ 2 - s_t .* s_a + s_a .^ 2) ; % Equivalent von mises stress 2D, [Pa]

YTS_margin = YTS ./ s_eqv - 1;

ep_a = (s_a - YTS) / E; % Plastic Strain Axial
ep_t = (s_t - YTS) / E; % Plastic Strain Tangential
ep_eqv = 2 / sqrt(3) .* sqrt(ep_t .^ 2 + ep_t .* ep_a + ep_a .^ 2); % Equivalent plastic strain from distortion energy theory, [MPa]

for i = 1:steps
    life_eqn = @(N) 3.5 * ft_ult / E * N ^ -0.12 + (e / N) ^ 0.6 - 2* (ep_eqv(i) + YTS / E);
    N(i) = fsolve(life_eqn, 5) / 4; 
end


%% Outputs

face_thermals = figure;
hold on 
% grid on 
set(gcf,'color','w')
set(gca,'TickLabelInterpreter','latex')

plot(htc_array, T_faceplate, 'Linewidth', 2)

dim = [0.5 0.4 0.2 0.2];
str = {
    'Wall Thickness [in]: ' + string(t_w / 0.0254),...
    'Fuel Velocity [m/s]: ' + string(round(velo,2)),...
    'Material Conductivity [W/m-K]: ' + string(k_w),...
 };
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex');

xlabel("Gas HTC [W/$m^2$-K]", 'Interpreter', 'latex')
ylabel("Wall Temp [K]", 'Interpreter', 'latex')
title("Injector Faceplate Steady-State Temperature",'Interpreter','latex')
exportgraphics(face_thermals, 'InjectorFaceThermals.png','Resolution',600)


face_stress = figure;
hold on 
% grid on 
set(gcf,'color','w')
set(gca,'TickLabelInterpreter','latex')

plot(T_faceplate, ep_eqv * 100, 'Linewidth', 2)

xlabel("Faceplate Temp [K]", 'Interpreter', 'latex')
ylabel("Equivalent Strain [\%]", 'Interpreter', 'latex')
title("Injector Faceplate Steady-State Plastic Strain",'Interpreter','latex')
