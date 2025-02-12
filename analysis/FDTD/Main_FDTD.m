clc
clear variables
% close all

[STC, modulation] = assemble_STC('Sinusoidal (discrete)', 3000);


%% Preliminary analysis of wave velocities

[x_grid, t_grid] = meshgrid(linspace(0, modulation.lambda, 100), linspace(0, modulation.period, 100));
[E_grid, J_grid, A_grid, rho_grid] = evaluate_structural_properties(x_grid, t_grid, STC);

c0_longitudinal = sqrt(E_grid ./ rho_grid);
c0_transverse   = sqrt((E_grid .* J_grid) ./ (rho_grid .* A_grid));

clear x_grid t_grid
clear E_grid J_grid A_grid rho_grid


%% Excitation settings

excitation = struct();
excitation.frequency  = 10.5 * 1e3;
excitation.amplitude  = 10 * 1e3;
excitation.coordinate = 50 / 100;

[~, final_time_force] = tone_burst(0, excitation.frequency, excitation.amplitude);


%% Numerical simulation settings

Nc = 50;

dx = STC(1).Beam.L - STC(1).Piezo.L;
dt = 0.2 * dx^2 / max(c0_transverse, [], 'all');
% dt = dx / max(c0_longitudinal, [], 'all');

x_max = Nc * modulation.lambda;
t_max = 3/9*7*(1 - excitation.coordinate) * x_max / (mean(c0_transverse, 'all') / dx);
% t_max = (1 - excitation.coordinate) * x_max / mean(c0_longitudinal, 'all');

x = 0:dx:x_max + dx * mod(ceil(x_max/dx), 2);
t = 0:dt:t_max + dt * mod(ceil(t_max/dt), 2);


%% Structural properties and force evaluation

[x_grid, t_grid] = meshgrid(x, t);
[E_grid, J_grid, A_grid, rho_grid] = evaluate_structural_properties(x_grid, t_grid, STC);

EJ_grid = E_grid .* J_grid;
rhoA_grid = rho_grid .* A_grid;

excitation.force = tone_burst(t, excitation.frequency, excitation.amplitude);


%% Finite difference finite time

tic
fprintf('FDFT simulation: \n');

try

    % U = solve_FDTD_rod(x, t, excitation, E_grid, rho_grid);
    U = solve_FDTD_beam(x, t, excitation, EJ_grid, rhoA_grid);

    if (any(isnan(U), 'all'))
        error('Simulation diverged.')
    end

catch ME
    fprintf(2, [ME.message '\n']);
    return
end

fprintf('\b%.2fs\n', toc);

% U_norm = norm( U(:, 1+round(excitation.coordinate * size(U, 2))) );


%% Fourier transform

Uz = zeros(size(U));
UU = U;
% UU = [Uz,Uz,Uz;Uz,U,Uz;Uz,Uz,Uz];

scale = size(UU, 1) / size(U, 1);
mu = (-scale*length(x)/2 : scale*length(x)/2-1) * modulation.lambda * 2*pi / (scale * max(x));
ff = (-scale*length(t)/2 : scale*length(t)/2-1) * 1 / (scale * max(t));

tic
fprintf('FFT decomposition: \n');
UU_fft = abs(fftshift(fft2(flipud(UU)))) * dt*dx;
fprintf('\b%.2fs\n', toc);


%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')

figure_fdtd = figure('Name', ['FDTD: ' modulation.label ' @' num2str(modulation.omega / (2*pi) * 1e-3, 3) 'kHz']);
tabgroup = uitabgroup(figure_fdtd);
colormap(jet(64))

% Create tabs
tab1 = uitab(tabgroup, 'Title', ['Structural properties | size([x,t]) = [' num2str(length(x)) ',' num2str(length(t)) ']']);
tab2 = uitab(tabgroup, 'Title', 'Force and Waterfall');
tab3 = uitab(tabgroup, 'Title', ['Dispersion diagram | [Nx,Nt,Nc]=[' num2str(length(x)) ',' num2str(length(t)) ',' num2str(Nc) ']']);


% Structural properties
axes('parent', tab1);
x_idxs = 1:min(round(modulation.lambda / diff(x(1:2))), length(x));
t_idxs = 1:min(round(modulation.period / diff(t(1:2))), length(t));
tile_EJ_grid = plot_structural_property(x_grid(t_idxs, x_idxs), t_grid(t_idxs, x_idxs), EJ_grid(t_idxs, x_idxs), 'Unit cell (EJ)', 'EJ [Nm^2]');
tile_rhoA_grid = plot_structural_property(x_grid(t_idxs, x_idxs), t_grid(t_idxs, x_idxs), rhoA_grid(t_idxs, x_idxs), 'Unit cell (rhoA)', 'rhoA [Kg/m]');
linkprop([tile_EJ_grid tile_rhoA_grid], {'View', 'XLim', 'YLim'});
clear x_idxs t_idxs

% Force and Waterfall
axes('parent', tab2);
% tile_force = plot_excitation_force(t, excitation.force);
tile_waterfall = plot_waterfall(x, t, U, 20);

% Dispersion diagram
axes('parent', tab3);
ff_idxs = find(ff >= 0 & ff <= 20e3);
mu_idxs = find(mu >= -5*pi & mu <= 5*pi);
plot_dispersion_diagram('fdtd', mu(mu_idxs), ff(ff_idxs), UU_fft(ff_idxs, mu_idxs), 10);
overlay(modulation.label, modulation.omega / (2*pi));
ylim([0 20])
grid minor

tabgroup.SelectedTab = tab3;