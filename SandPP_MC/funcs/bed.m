function [vz, vr, eps, CF] = bed(T, FI, z)
% Copyright 2024, All Rights Reserved
% Code by Stefano Iannello
% For Paper, "The behaviour of plastic particles during pyrolysis in 
%        bubbling fluidized bed reactors: Incipient agglomeration and 
%        axial segregation"
% by S. Iannello, A. Sebastiani, M. Errigo, M. Materazzi


T_amb = 25 + 273.15;                                                        % Ambient temperature [K]
g = 9.81;                                                                   % Acceleration of gravity [m/s^2]


%% Emulsion phase
if T == 500 + 273.15
    Qmf_amb = 13.4;                                                         % Experimental minimum fluidization flow rate at normal conditions [Nlpm]
elseif T == 600 + 273.15
    Qmf_amb = 10.4;
elseif T == 650 + 273.15
    Qmf_amb = 10;
else
    error('Temperature should be 500, 600 or 650 degC')
end

d_bed = 0.14;                                                               % Bed diameter [m]
A_bed = pi * d_bed^2/4;                                                     % Bed cross sectional area [m^2]
umf = Qmf_amb / (A_bed * 1000 * 60) * T/T_amb;                              % Recalculate the minimum fluidization velocity at operating conditions using the ideal gas law [m/s]
u = FI * umf;                                                               % Fluidization velocity [m/s]
ds = 250 * 10^-6;                                                           % Mean diameter of bed material [m]

if T == 500 + 273.15
    eps_mf = 0.4649;                                                        % Bed voidage at experimental minimum fluidization conditions
     if FI == 1.0
        rho_bed = 1418;                                                     % Bed bulk density at experimental condition, [kg/m^3]
        eps = eps_mf;                                                       % Bed voidage at experimental conditions
     elseif FI == 1.25
        rho_bed = 1331;
        eps = 0.4977;
     elseif FI == 2.0
        rho_bed = 1280;
        eps = 0.5170;
     end
elseif T == 600 + 273.15
    eps_mf = 0.4755;
    if FI == 1.0
        rho_bed = 1390;
        eps = eps_mf; 
    elseif FI == 1.25
        rho_bed = 1334;
        eps = 0.4966;
    elseif FI == 2.0
        rho_bed = 1287;
        eps = 0.5143;
    end
elseif T == 650 + 273.15
    eps_mf = 0.4630;
    if FI == 1.0
        rho_bed = 1423;
        eps = eps_mf;
    elseif FI == 1.25
        rho_bed = 1352;  
        eps = 0.4898;
    elseif FI == 2.0
        rho_bed = 1307;
        eps = 0.5068;
    end
end

rhos = 2650;                                                                % Particle density of bed solids [kg/m^3]
rho_bulk = 1650;                                                            % Bulk density of bed at ambient conditions [kg/m^3]
eps_fix = (rhos - rho_bulk) / rhos;                                         % Voidage of fixed bed
ue = umf / eps_mf;                                                          % Emulsion gas velocity [m/s] (from from Kunii and Levenspiel, 1991)


%% Bubble phase (correlations from from Kunii and Levenspiel, 1991, unless otherwise specified)
fw = 0.25;                                                                  % Bubble wake fraction for Geldart B 
db0 = (2.78 / g) * (u - umf)^2;                                             % Initial bubble diameter [m]
dbm = (0.65 * ((pi/4) * (d_bed * 100)^2 * (u*100 - umf*100))^0.4) / 100;    % Maximum bubble diameter [m] (multipl and division by 100 as the correlation gives reult in cm)
db = dbm - (dbm - db0) * exp(-0.3 * z/d_bed);                               % Bubble diameter along bed height [m]
ubr = 0.711 * (g * db).^0.5;                                                % Rise velocity of a single bubble [m/s]
ub = (u - umf) + ubr;                                                       % Bubbles rise velocity along bed height [m/s]
us_up = ub;                                                                 % Upward velocity of wake solids [m/s]

if ub < ue                                                                  % Slow bubbles
    delta = (u - umf) / (ub + 2*umf);
elseif ue < ub && ub < 5*ue                                                 % Intermediate bubbles with thick clouds
    delta = ((u - umf) / (ub + umf) + (u - umf) / ub) / 2;
elseif ub > 5 * ue                                                          % Fast bubbles
    delta = (u - umf) / (ub - umf);
else
    delta = u / ub;                                                         % Vigorous bubbling
end


%% Bed solids motion (correlations from from Kunii and Levenspiel, 1991, unless otherwise specified)
us_down = fw * delta * ub / (1 - delta - fw*delta);                         % Velocity of sinking bed solids [m/s]
vz = us_up - us_down;                                                       % Net average vertical velocity of bed solids [m/s]
alpha = (1 + 0.77) / 2;                                                     % Parameter depending on Geldart classification. For A/AB is 1, for BD is 0.77.
D_sh = 3/16 * (delta / (1 - delta)) * alpha^2 * db * ubr * (((ubr + 2*ue) ...
    / (ubr - ue))^(1/3) - 1);                                               % Horizontal dispersion coefficient of bed particles [m^2/s]
vr = 3 * sqrt(8) * D_sh * (1 - eps) / ds;                                   % Average horizontal velocity of bed solids [m/s] (from Migliozzi et al., 2017) 
CF = 56400*((1 - eps)/(1 - eps_fix))^2 * (1 - (1 - eps)/(1 - eps_fix)) * u; % Collision frequency of bed solids [Hz] (from Buffiere and Moletta, 2000)
end