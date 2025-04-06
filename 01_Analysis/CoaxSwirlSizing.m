%% Coaxial Swirl Injector Sizer
% % Author: Alex Suppiah
% % First Created: 2/12/2025
% V. Bazarov, V. Yang, P. Puri
% Liquid Atomization

close all
clear
clc

%% Engine Parameters
fuel = 'Ethanol'; % fuel definition
oxidizer = 'O2'; % oxidizer definition
Tl_fu = 293.15; % [K]
Tl_ox = 90.17; % [K]
Pc = 250 * 6894.76; % Chamber pressure, [Pa]
OF = 1.2;
mdot_ox = 0.673131077; % Massflow ox, [kg/s]
mdot_fu = 0.560640169; % Massflow fuel, [kg/s]
T_tot = 550; % Total thrust, [lbf]


%% Injector Parameters
% Design parameters
element_count = 6;

Cd_orifices = 0.7; % Circular orifice Cd (machined 0.7, printed 0.7)
dP_pct_ox = 0.2; % Ox stiffness
dP_pct_fu = 0.2; % Fuel stiffness

% Calculations
Pman_ox = Pc * (1 + dP_pct_ox); % Ox manifold pressure
Pman_fu = Pc * (1 + dP_pct_fu); % Fuel manifold pressure

dP_ox = Pman_ox - Pc; % Ox side dP
dP_fu = Pman_fu - Pc; % Fuel side dP

mdot_ox_elem = mdot_ox / element_count; % Ox mdot per element, [kg/s]
mdot_fu_elem = mdot_fu / element_count; % Fu mdot per element, [kg/s]
T_per_elem = T_tot / element_count; % Thrust per element, [lbf]


%% Element Geometric Parameters
RN = 1.25; % Recess number 
total_length = 0.06; % Total element length


%% Inner Swirler Parameters
i_ox = 4; % Ox side tangential swirl inlets
alpha_inner = 75; % Desired spray angle, [deg]
ds_do = 1; % Ratio of swirl arm to outlet radius, inner (Recommend 2-5)
ds_do_outer = 0.6; % Ratio of swirl arm to outlet radius, outer (Recommend 2-5)
t_wall = 2 / 1000; % Ox post wall thickness, [m]


%% Outer Swirler Parameters
% Design Parameters
t_gap = 0.3 / 1000; % Gas core size uggested value from Bazarov 0.3mm, [m]

% Geometry Parameters
i_fu = 4; % Fuel side tangential swirl inlets
t_wall_outer = 3 / 1000; % Fuel wall thickness (unimportant parameter), [m]


%% Propellant Properties
rho_ox = py.CoolProp.CoolProp.PropsSI('D','T', Tl_ox,'P', Pman_ox, oxidizer); % Ox density in manifold, [kg/m^3]
rho_fu = 785; %py.CoolProp.CoolProp.PropsSI('D','T', Tl_fu,'P', Pman_fu, fuel); % Fu density in manifold, [kg/m^3]
mu_ox = py.CoolProp.CoolProp.PropsSI('V','T', Tl_ox,'P', Pman_ox, oxidizer); % Ox dynamic viscosity, [Pa-s]
mu_fu = py.CoolProp.CoolProp.PropsSI('V','T', Tl_fu,'P', Pman_fu, fuel); % Fu dynamic viscosity, [Pa-s]


%% Inner Swirler Sizing - Ox
options = optimoptions('fsolve','Display','none');

% Guess filling efficiency K value to find e (coefficient of passage fullness)
K_guess = 2; 
error_inner = alpha_inner; 
tolerance_inner = 0.01;
i = 0;

while abs(error_inner) > tolerance_inner && i < 100
    % Solve for e implicitly
    K_func = @(e) (1 - e) * sqrt(2) / (e * sqrt(e)) - K_guess; % Eq. 5-65, Beyvel and Orzechowski 
    e = fzero(K_func, 1); % Coefficient of passage fullness
    
    % Find Cd
    Cd = e * sqrt(abs(e / (2 - e))); % Discharge Coefficient (Eq. 5-66, Beyvel and Orzechowski)
        
    % Find Outlet Diameter
    d_o = sqrt(4 * mdot_ox_elem / (pi * Cd * sqrt(2 * rho_ox * dP_ox))); % Outlet Diameter (Eq. 5-82, Beyvel and Orzechowski) 
    
    % Find Inlet Diameter
    R = ds_do * d_o / 2;
    d_i = sqrt(2 * R * d_o / (i_ox * K_guess)); % Inlet diameter (Eq. 5-84, Beyvel and Orzechowski)
    
    % Find Re
    Re = 4 * mdot_ox_elem / (pi * mu_ox * sqrt(i_ox) * d_i); % Reynolds number
    
    % Find empirical friction factor, frictional loss coeff in inlets
    lambda = exp(25.8 / log(Re) ^ 2.58 - 2); % (Eq. 5-87, Beyvel and Orzechowski)
    
    % Find new K accounting for frictional losses in inlets
    K_lambda = R * (d_o / 2) / (i_ox * (d_i / 2) ^ 2 + lambda / 2 * R * (R - d_o / 2));  % (Eq. 5-81, Beyvel and Orzechowski)
    
%         % Find new filling efficiency with K_lambda
%         K_func = @(e) (1 - e) * sqrt(2) / (e * sqrt(e)) - K_lambda; % Eq. 5-65, Beyvel and Orzechowski 
%         e = fzero(K_func, 0.7); % Coefficient of passage fullness
%         
%         % Find Cd
%         Cd = e * sqrt(abs(e / (2 - e))); % Discharge Coefficient (Eq. 5-66, Beyvel and Orzechowski)
% 
%         % Find new Cd taking momentum loss into account
%         Cd = Cd / sqrt(1 + Cd ^ 2 * K_lambda ^ 2 / R ^ 2); % Discharge Coefficient (Eq. 5-66, Beyvel and Orzechowski)
% 
%         % Find new Outlet Diameter
%         d_o = sqrt(4 * mdot_ox_elem / (pi * Cd * sqrt(2 * rho_ox * dP_ox))); % Outlet Diameter (Eq. 5-82, Beyvel and Orzechowski) 
% 
%         % Find new R
%         R = ds_do * d_o / 2;
% 
%         % Find new K_guess
%         K_guess = 2 * R * d_o / (i_ox * d_i ^ 2);

    % Solve for S implicitly
    Cd_func = @(S) sqrt(abs(1 - Cd ^ 2 * K_guess ^ 2)) - S * sqrt(abs(S ^ 2 - Cd ^ 2 * K_guess ^ 2)) ...
        - Cd ^ 2 * K_guess ^ 2 * log((1 + sqrt(abs(1 - Cd ^ 2 * K_guess ^ 2))) / (S + sqrt(abs(S ^ 2 - Cd ^ 2 * K_guess ^ 2)))) - Cd; 
    S = fzero(Cd_func, 1); % Ratio gas core diameter to swirl outlet diameter (Eq. 5-58, Beyvel and Orzechowski)

    % Find viscous cone angle
    alpha_actual = atand(2 * Cd * K_lambda / sqrt(abs((1 + S) ^ 2 - 4 * Cd ^ 2 * K_guess ^ 2))) * 2; % (Eq. 5-75, Beyvel and Orzechowski)
    
    K_guess = K_lambda * alpha_inner / alpha_actual;

    error_inner = alpha_inner - alpha_actual;
    
    i = i + 1; 
end

if error_inner > tolerance_inner
    fprintf("Inner swirler did not converge, reduce target alpha. Error: %0.2f\n", error_inner)
end

area_i_design = pi / 4 * d_i ^ 2;
area_i = area_i_design / Cd_orifices; 

% Inner Swirler Parameters
d_o_inner = d_o; % [m]
d_i_inner = sqrt(area_i * 4 / pi); % [m] 
Cd_inner = Cd; 
K_inner = K_guess; 
alpha_inner = alpha_actual;
swirl_arm = ds_do * d_o_inner; % [m]


%% Outer Swirler Sizing - Fu
t_f_guess = 0.001; 
error_outer = t_f_guess;
tolerance_outer = 0.00001; 

i = 0;
while abs(error_outer) > tolerance_outer && i < 100
    % Find outer swirl outlet diam
    d_o = d_o_inner + 2 * t_gap + 2 * t_wall + 2 * t_f_guess; 
    
    % Find flow areas
    A_eff = mdot_fu / sqrt(2 * rho_fu * dP_fu);
    A_o = pi / 4 * d_o ^ 2; 
    
    % Find Cd 
    Cd = A_eff / A_o; 
    
    % Solve for e implicitly
    Cd_func = @(e) e * sqrt(abs(e / (2 - e))) - Cd; % Discharge Coefficient (Eq. 5-66, Beyvel and Orzechowski)
    e = fzero(Cd_func, 1); % Coefficient of passage fullness
    
    % Find K
    K = (1 - e) * sqrt(2) / (e * sqrt(e)); % Eq. 5-65, Beyvel and Orzechowski 
    
    % Solve for S implicitly
    Cd_func = @(S) sqrt(abs(1 - Cd ^ 2 * K ^ 2)) - S * sqrt(abs(S ^ 2 - Cd ^ 2 * K ^ 2)) ...
        - Cd ^ 2 * K ^ 2 * log((1 + sqrt(abs(1 - Cd ^ 2 * K ^ 2))) / (S + sqrt(abs(S ^ 2 - Cd ^ 2 * K ^ 2)))) - Cd; 
    S = fzero(Cd_func, 1); % Ratio gas core diameter to swirl outlet diameter (Eq. 5-58, Beyvel and Orzechowski)
    
    % Find Inlet Diameter
    R = ds_do_outer * d_o / 2;
    d_i = sqrt(2 * R * d_o / (i_fu * K)); % Inlet diameter (Eq. 5-84, Beyvel and Orzechowski)
    
    % Find Re
    Re = 4 * mdot_fu_elem / (pi * mu_fu * sqrt(i_fu) * d_i); % Reynolds number
    
    % Find empirical friction factor, frictional loss coeff in inlets
    lambda = exp(25.8 / log(Re) ^ 2.58 - 2); % (Eq. 5-87, Beyvel and Orzechowski)
    
    % Find new K accounting for frictional losses in inlets
    K_lambda = R * (d_o / 2) / (i_fu * (d_i / 2) ^ 2 + lambda / 2 * R * (R - d_o / 2));  % (Eq. 5-81, Beyvel and Orzechowski)

    % Find viscous cone angle
    alpha_actual = atand(2 * Cd * K_lambda / sqrt(abs((1 + S) ^ 2 - 4 * Cd ^ 2 * K ^ 2))) * 2; % (Eq. 5-75, Beyvel and Orzechowski)
    
    % Find film thickness
    t_f_actual = 1 / 2 * (d_o - S * d_o); 
    
    error_outer = t_f_guess - t_f_actual;
    t_f_guess = t_f_actual; 

    i = i + 1; 
end

if error_outer > tolerance_outer
    fprintf("Outer swirler did not converge. \n")
end

area_i_design = pi / 4 * d_i ^ 2;
area_i = area_i_design / Cd_orifices; 

% Outer Swirler Parameters
d_o_outer = d_o; % [m]
d_i_outer = sqrt(area_i * 4 / pi);
Cd_outer = Cd; 
K_outer = K; 
t_f = t_f_actual;
alpha_outer = alpha_actual; 
swirl_arm_outer = ds_do_outer * d_o_outer; % [m]


%% Total Swirler Parameters
velo_tot_inner = mdot_ox_elem / (pi / 4 * d_o_inner ^ 2 * rho_ox); % [m/s]
velo_axial_inner = velo_tot_inner * cosd(alpha_inner / 2); % [m/s]
velo_tan_inner = velo_tot_inner * sind(alpha_inner / 2); % [m/s]

velo_tot_outer = mdot_fu_elem / (pi / 4 * d_o_outer ^ 2 * rho_fu); % [m/s]
velo_axial_outer = velo_tot_outer * cosd(alpha_outer / 2); % [m/s]
velo_tan_outer = velo_tot_outer * sind(alpha_outer / 2); % [m/s]

MR_axial = mdot_ox_elem * velo_axial_inner + mdot_fu_elem * velo_axial_outer; % [kg-m/s]
MR_tan = mdot_ox_elem * velo_tan_inner + mdot_fu_elem * velo_tan_outer; % [kg-m/s]
MR_tot = sqrt(MR_axial ^ 2 + MR_tan ^ 2); % [kg-m/s]
alpha_tot = atand(MR_tan / MR_axial) * 2; % [deg]

Lc = (d_o_outer / 2 - d_o_inner / 2) / tand(alpha_inner / 2); % Spray cone axial length, [m]
Lr = Lc * RN; % Recess length, [m]
gap = d_o_outer / 2 - d_o_inner / 2 - t_wall; % Space between posts, [m]
recess_length = Lc; 


%% Display Results
n = 2;
unit_convert = 1000; % m to mm 
plot_length = (total_length + total_length/5) * unit_convert;

% Inner swirler
ox_length = total_length - Lr;
r_o_inner = d_o_inner / 2;
r_i_inner = d_i_inner / 2; 
inlet_line_i = linspace(0, r_o_inner, n); 
inner_wall_i = linspace(r_o_inner, r_o_inner, n); 
outer_wall_i = linspace(r_o_inner + t_wall, r_o_inner + t_wall, n); 
post_thickness = linspace(r_o_inner, r_o_inner + t_wall, n); 
spray_cone_i = linspace(r_o_inner, r_o_inner + tand(alpha_inner / 2) * Lr, n);

% Outer swirler
fu_start = 0.02; 
fu_length = total_length;
r_o_outer = d_o_outer / 2;
r_i_outer = d_i_outer / 2; 
inlet_line_o = linspace(r_o_inner + t_wall, r_o_outer, n); 
inner_wall_o = linspace(r_o_outer, r_o_outer, n); 
outer_thickness = linspace(r_o_outer, r_o_outer + t_wall_outer, n); 
outer_wall_o = linspace(r_o_outer + t_wall_outer, r_o_outer + t_wall_outer, n); 
end_wall = linspace(r_o_outer, r_o_outer + t_wall_outer, n); 
spray_cone_o = linspace(r_o_outer, r_o_outer + tand(alpha_outer / 2) * Lr, n);
spray_cone_t = linspace(r_o_outer, r_o_outer + tand(alpha_tot / 2) * total_length/5, n);

% Plot
element_plot = figure(1);
hold on 
% grid on 
axis equal
set(gcf,'color','w')
set(gca,'TickLabelInterpreter','latex')

% Ox post
plot(linspace(0, 0, n), inlet_line_i * unit_convert, 'k', 'Linewidth', 2) % Plot inlet surface
plot(linspace(0, ox_length, n) * unit_convert, inner_wall_i * unit_convert, 'k', 'Linewidth', 2) % Plot inner wall
plot(linspace(0, ox_length, n) * unit_convert, outer_wall_i * unit_convert, 'k', 'Linewidth', 2) % Plot outer wall
plot(linspace(ox_length, ox_length, n) * unit_convert, post_thickness * unit_convert, 'k', 'Linewidth', 2) % Plot outer wall

plot(linspace(0, 0, n), -inlet_line_i * unit_convert, 'k', 'Linewidth', 2) % Plot inlet surface
plot(linspace(0, ox_length, n) * unit_convert, -inner_wall_i * unit_convert, 'k', 'Linewidth', 2) % Plot inner wall
plot(linspace(0, ox_length, n) * unit_convert, -outer_wall_i * unit_convert, 'k', 'Linewidth', 2) % Plot outer wall
plot(linspace(ox_length, ox_length, n) * unit_convert, -post_thickness * unit_convert, 'k', 'Linewidth', 2) % Plot outer wall

% Fuel post
plot(linspace(fu_start, fu_start, n) * unit_convert, inlet_line_o * unit_convert, 'r', 'Linewidth', 2) % Plot inlet surface
plot(linspace(fu_start, fu_length, n) * unit_convert, inner_wall_o * unit_convert, 'r', 'Linewidth', 2) % Plot inner wall
plot(linspace(fu_start, fu_length, n) * unit_convert, outer_wall_o * unit_convert, 'r', 'Linewidth', 2) % Plot outer wall
plot(linspace(fu_length, fu_length, n) * unit_convert, end_wall * unit_convert, 'r', 'Linewidth', 2) % Plot outer wall

plot(linspace(fu_start, fu_start, n) * unit_convert, -inlet_line_o * unit_convert, 'r', 'Linewidth', 2) % Plot inlet surface
plot(linspace(fu_start, fu_length, n) * unit_convert, -inner_wall_o * unit_convert, 'r', 'Linewidth', 2) % Plot inner wall
plot(linspace(fu_start, fu_length, n) * unit_convert, -outer_wall_o * unit_convert, 'r', 'Linewidth', 2) % Plot outer wall
plot(linspace(fu_length, fu_length, n) * unit_convert, -end_wall * unit_convert, 'r', 'Linewidth', 2) % Plot outer wall

% Spray cones
plot(linspace(ox_length, fu_length, n) * unit_convert, spray_cone_i * unit_convert, 'b', 'Linewidth', 1.5)
plot(linspace(fu_length, fu_length + Lr, n) * unit_convert, spray_cone_o * unit_convert, 'Linewidth', 1.5, 'Color', [0.9290 0.6940 0.1250])
plot(linspace(fu_length * unit_convert, plot_length, n), spray_cone_t * unit_convert, 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980])

plot(linspace(ox_length, fu_length, n) * unit_convert, -spray_cone_i * unit_convert, 'b', 'Linewidth', 1.5)
plot(linspace(fu_length, fu_length + Lr, n) * unit_convert, -spray_cone_o * unit_convert, 'Linewidth', 1.5, 'Color', [0.9290 0.6940 0.1250])
plot(linspace(fu_length * unit_convert, plot_length, n), -spray_cone_t * unit_convert, 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980])

% Orifices
theta = linspace(0,2*pi,100);
plot((r_i_inner + r_i_inner*cos(theta)) * unit_convert, (swirl_arm + r_i_inner*sin(theta)) * unit_convert, 'Linewidth', 1.5) % Ox orifice
plot((r_i_outer + fu_start + r_i_outer*cos(theta)) * unit_convert, (swirl_arm_outer + r_i_outer*sin(theta)) * unit_convert, 'Linewidth', 1.5) % Fuel orifice

% Annotations
dim = [0.2 0.15 0.2 0.2];
str = {
    'Inner Swirler (' + string(i_ox) + ' inlets)',...
    'Orifice Diam [mm]: ' + string(round(d_i_inner * unit_convert,2)),...
    'Outlet Diam [mm]: ' + string(round(d_o_inner * unit_convert,2)),...
    'Swirl Arm [mm]: ' + string(round(swirl_arm * unit_convert,2)),...
    'Post Thickness [mm]: ' + string(round(t_wall * unit_convert,2)),...
 };
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex');

dim = [0.55 0.15 0.2 0.2];
str = {
    'Outer Swirler (' + string(i_fu) + ' inlets)',...
    'Orifice Diam [mm]: ' + string(round(d_i_outer * unit_convert,2)),...
    'Outlet Diam [mm]: ' + string(round(d_o_outer * unit_convert,2)),...    
    'Annular Gap [mm]: ' + string(round(gap * unit_convert,2)),...
    'Recess Length [mm]: ' + string(round(recess_length * unit_convert,2)),...
 };
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex');

dim = [0.35 0.72 0.2 0.2];
str = {
    'Element Count: ' + string(element_count),...
    'Thrust/Element [lbf]: ' + string(round(T_per_elem)),...
    'ds/do: ' + string(ds_do),...
    'Total Length [mm]: ' + string(total_length * unit_convert),...
    'Inner Cone Angle [deg]: ' + string(round(alpha_inner,1)),...
    'Outer Gas Core Size [mm]: ' + string(t_gap * unit_convert)
 };
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex');

xlabel("Length [mm]", 'Interpreter', 'latex')
ylabel("Width [mm]", 'Interpreter', 'latex')
title("Liquid-Liquid Coaxial Bi-Swirl Element",'Interpreter','latex')
xlim([0, plot_length])
exportgraphics(element_plot,'CoaxSwirlSizing.png','Resolution',600)

