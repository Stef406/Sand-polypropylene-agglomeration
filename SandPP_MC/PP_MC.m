% Copyright 2024, All Rights Reserved
% Code by Stefano Iannello
% For Paper, "The behaviour of plastic particles during pyrolysis in 
%        bubbling fluidized bed reactors: Incipient agglomeration and 
%        axial segregation"
% by S. Iannello, A. Sebastiani, M. Errigo, M. Materazzi

clc; clear; close all;
addpath('./funcs');


%% Simulation parameters
T = input('T Bed (Celsius) = ') + 273.15;                                   % Bed temperature [K]
FI = input('U/Umf = ');                                                     % Fluidization index [U/Umf]
x_O2 = input('O2 concentration = ');                                        % Molar fraction of oxygen in reaction environment (0 if pyrolysis, 0.21 if combustion in air)
z = input('Height of the PP particle in bed (cm) = ') / 100;                % Distance of the PP particle from the distributor plate within the fluidized bed reactor [m]
d_pp0 = 12 * 10^-3;                                                         % Initial diameter of polypropylene particle [m]
[vz, vr, eps, CF] = bed(T, FI, z);                                    
V_pp0 = pi / 6 * d_pp0^3;                                                   % Initial volume of PP particle [m^3]
l = d_pp0;                                                                  % Size simulation box (cube) [m]
V_box = l^3;                                                                % Volume of the simulation box [m^3]
V_sinbox = (V_box - V_pp0) * (1 - eps);                                     % Intitial Volume occupied by the sand in the simulation box [m^3]
d_bm = 250 * 10^-6;                                                         % Mean diameter of bed material particles [m]
rho_bm = 2650;                                                              % Density of bed material particles [kg/m^3]
V_s1 = pi / 6 * d_bm^3;                                                     % Volume of a single sand particle [m^3]
M_sinbox = rho_bm * V_sinbox;                                               % Mass of sand in the simulation box [g]
n_particlesMax = V_sinbox / V_s1;                                           % Maximum number of sand particles allowed in the simulation box
n_particles = round(n_particlesMax);                                        % Number of particles in the simulation box
m_pp0 = 0.65 * 10^-3;                                                       % Initial mass of polypropylene particle [kg]
rho_pp = 697;                                                               % Density of polypropylene [kg/m^3]
t_sim = 0.1;                                                                % Simulation time [s]
dt = 1/CF;                                                                  % Time step [s]
time = 0:dt:t_sim;                                                          % Time vector [s]
k = devol(T, x_O2, d_pp0);                                                  % Devolatilization rate constant [1/s]


%% Initialize random location and velocities of sand particles
rng(92)                                                                     % Set random number generator with seed for reproducibility

pos = -l/2 + (l/2 + l/2) * rand(n_particles, 3);                            % Initialize position of bed material particles
vel = [vr, vr, vz] .* ones(n_particles, 3);                                 % Velocity of bed material particles
pos_pp = [l/2 l/2 l/2];                                                     % Position of polypropylene particle = centre of the simulation box
prog_bar = waitbar(0);                                                      % Set up progress bar


%% Initialize video
video = VideoWriter('Agglomeration');                                       % Open video file
video.FrameRate = 6;                                                        % Adjust as you like
open(video)


%% Simulation loop
n_agglomerated = 0;                                                         % Counter for number of particles stuck to polypropylene
n_aggl = zeros(1, n_particles, length(time));                               % Store number of sand particles stuck on the PP
n_agglpcycle = zeros(1, length(time));                                      % Store number of sand particles stuck on the PP per simulation cycle
dist = zeros(1, n_particles, length(time));                                 % Store distance between sand particle and centre of PP [m]
dist_aggl = zeros(1, n_particles, length(time));                            % Store distance between stuck sand particle and centre of PP [m]
d_pp = zeros(1, length(time));                                              % Store diameter of shrinking of PP [m]
r_pp = zeros(1, length(time));                                              % Store radius of shrinking of PP [m]
cycle = 1;                                                                  % Counter number of simulation cycle



tic
for t = 0:dt:t_sim

    % Shrink polypropylene particle (pseudo-first order kinetic)
    d_pp(cycle) = d_pp0 * exp(-k * t / 3);
    r_pp(cycle) = d_pp(cycle) / 2;

    % Check for collisions between particles
    for i = 1 : n_particles

        % Calculate distance between particles
        dist(:, i, cycle) = sqrt(sum((pos(i, :) - pos_pp).^2));

            % Check if the bed material particle is within the shrinking PP sphere
            if  dist(:, i, cycle) <= r_pp(cycle) + d_bm/2 && dist(:, i, cycle) >= r_pp(cycle) - d_bm/2

                % Stick bed material particle to polypropylene
                n_agglomerated = n_agglomerated + 1;
                n_aggl(:, i, cycle) = n_agglomerated;
                n_agglpcycle(cycle) = length(find(n_aggl(:, :, cycle) > 0));
                pos_agglomerated(n_aggl(:, i, cycle), :) = pos(i, :);
                dist_aggl(:, i, cycle) = sqrt(sum((pos_agglomerated(n_aggl(:, i, cycle), :) - pos_pp).^2));
            end
    end


    % Update positions of bed material particles
    pos = pos + vel*dt;
    
   
    for i = 1:n_particles
        for j = 1:3

            % Reflect particles at walls of simulation box
            if pos(i, j) - d_bm/2 < 0
                pos(i, j) = pos(i, j) + l;
            elseif pos(i, j) + d_bm/2 > l
                pos(i, j) = pos(i, j) - l;
            end
        end
    end


    % Update progress bar
    waitbar(t/t_sim, prog_bar, sprintf('Simulating agglomeration... %.1f%%', 100 * t/t_sim));
    
    % Generate 3D graph of agglomeration in real-time (this visulization can affect performance. Non reccommended for long simulations)
    if t_sim <= 2
  
        figure(1)
        clf

        % Bed material particles
        scatter3(pos(:, 1) / max(pos(:, 1)), pos(:, 2) / max(pos(:, 2)), pos(:, 3) / max(pos(:, 3)), 3, 'filled','MarkerFaceColor', [0.86, 0.78, 0.63])
        hold on
        scatter3(pos_agglomerated(:, 1)/max(pos_agglomerated(:, 1)), pos_agglomerated(:, 2)/max(pos_agglomerated(:, 2)), pos_agglomerated(:, 3)/max(pos_agglomerated(:,3)), 3, 'filled','MarkerFaceColor', [0.4, 0.3, 0.1])
        
        % Polypropylene particle
        [x2, y2, z2] = sphere(100);
        x2 = x2 * r_pp(cycle)/l + pos_pp(1)/l;
        y2 = y2 * r_pp(cycle)/l+ pos_pp(2)/l;
        z2 = z2 * r_pp(cycle)/l + pos_pp(3)/l;
        surf(x2, y2, z2, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', [0, 0, 0.5])
        grid off
        axis equal
        xlim([0 l/l])
        ylim([0 l/l])
        zlim([0 l/l])
        xlabel('X [-]')
        ylabel('Y [-]')
        zlabel('Z [-]')
        title(sprintf('%.2f seconds', t))

        drawnow
    else

    end

    % Update simulation cycle number
    cycle = cycle + 1;

    % Record video frame
    frame = getframe(gcf);
    writeVideo(video, frame);
end
toc

% Close the video
close(video)

% Close progress bar
close(prog_bar)


%% Calculate the density of the agglomerate
cycle_p = 1:cycle - 1;

% without lead tracer
M_agg = zeros(1, length(n_agglpcycle));
V_agg = zeros(1, length(n_agglpcycle));
rho_agg = zeros(1, length(n_agglpcycle));

for i = 1 : length(n_agglpcycle)
    M_agg(i) = pi/6 * rho_pp * d_pp(i).^3 + n_agglpcycle(i) * rho_bm * V_s1;
    V_agg(i) = pi/6 * d_pp(i).^3 + n_agglpcycle(i) * V_s1;
    rho_agg(i) =  M_agg(i) ./ V_agg(i);
end

% with lead tracer
M_agg_l = zeros(1, length(n_agglpcycle));
V_agg_l = zeros(1, length(n_agglpcycle));
rho_agg_l = zeros(1, length(n_agglpcycle));
ml = 0.17 * 10^-3;                                                          % Mass of lead tracer within the PP particle [kg]

for i = 1 : length(n_agglpcycle)
    M_agg_l(i) = pi/6 * rho_pp * d_pp(i).^3 + n_agglpcycle(i) * rho_bm * V_s1 + ml;
    V_agg_l(i) = pi/6 * d_pp(i).^3 + n_agglpcycle(i) * V_s1;
    rho_agg_l(i) =  M_agg_l(i) ./ V_agg_l(i);
end

% PP density considering no agglomeration
m_pp = m_pp0 .* exp(-k .* time);
V_pp = pi/6 * d_pp.^3;
rho_pp_l = (m_pp + ml) ./ V_pp;

fprintf('Number of particles stuck to polypropylene: %d\n', n_agglomerated);
fprintf('Final density of PP-sand agglomerate: %d\n', rho_agg(end));

%% Plots
figure(2)
yyaxis left
semilogx(cycle_p, n_agglpcycle, 'k-', 'LineWidth', 1)
ylabel('N_{agg}')
hold on
yyaxis right
semilogx(cycle_p, d_pp * 1000, 'r-', 'LineWidth', 1)
ylabel('d_{p} [mm]')
xlabel('Number of cycles')
title('Sand particles agglomerated and diameter of shrinking PP', 'FontSize', 9)
legend({'N_{agg}','d_{p}'}, 'Location', 'south')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
hold off

figure(3)
plot(time, rho_agg, 'k-', 'LineWidth', 1)
xlabel('t [s]')
ylabel('\rho_{agg} [kg/m^{3}]')
title('Density evolution of the sand-PP agglomerate without lead tracer', 'FontSize', 9)

figure(4)
plot(time, rho_agg_l, 'k-', 'LineWidth', 1)
hold on 
plot(time, rho_pp_l, 'k--', 'LineWidth', 1)
xlabel('t [s]')
ylabel('\rho [kg/m^{3}]')
title('Density evolution of the sand-PP agglomerate and PP particle with lead tracer',  'FontSize', 9)
legend({'MC simulation (agglomerate density)','PP particle density'}, 'Location', 'best')
hold off