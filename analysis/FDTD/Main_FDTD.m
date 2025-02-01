% clc
clear variables
% close all

[STC, modulation] = assemble_STC('OFF-OFF-OFF');

E_mean = mean(arrayfun(@(unit) mean(unit.E{1}(linspace(1, modulation.period, 100))), STC));
rho_mean = mean(arrayfun(@(unit) mean(unit.rho{1}), STC));

% We suppose slightly higher velocity to avoid boundaries wave reflection
c0 = 1.1 * sqrt(E_mean / rho_mean);

%% Excitation settings
% The smaller fe or dfe, the longer the time needed

excitation = struct();
excitation.frequency  = 35e+03; %8.4 * 1e3;
excitation.amplitude  = 5e+03; %8.0 * 1e3;
excitation.coordinate = 50 / 100;


%% Numerical simulation settings

Nc = 100;

x_max = Nc * modulation.lambda;
t_max = (1 - excitation.coordinate) * x_max / c0;

Nx = floor(12 * 3 * Nc);
% Nt = floor(sqrt(2) * Nx);
Nt = 1.5 * Nx * (t_max / x_max * c0);

x = linspace(0, x_max, Nx - mod(Nx, 2));
t = linspace(0, t_max, Nt - mod(Nt, 2));


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
    U = solve_FDTD_beam(x, t, excitation, E_grid, rho_grid, J_grid, A_grid);
catch ME
    disp(ME.message);
    return
end

fprintf('\b%.2fs\n', toc);

U_norm = norm( U(:, 1+round(excitation.coordinate * size(U, 2))) );


%% Fourier transform

Uz = zeros(size(U));
UU = U;
% UU = [Uz,Uz,Uz;Uz,U,Uz;Uz,Uz,Uz];

scale = size(UU, 1) / size(U, 1);

mu = (-scale*length(x)/2 : scale*length(x)/2-1) * modulation.lambda * 2*pi / (scale * max(x));
ff = (-scale*length(t)/2 : scale*length(t)/2-1) * 1 / (scale * max(t));

tic
fprintf('FFT decomposition: \n');
UU_fft = abs(fftshift(fft2(UU)));
fprintf('\b%.2fs\n', toc);


%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')

figure_fdtd = figure('Name', ['FDTD: ' modulation.label ' @' num2str(modulation.omega / (2*pi) * 1e-3, 3) 'kHz']);
tabgroup = uitabgroup(figure_fdtd);
colormap(jet(64))


% % Structural properties
% axes('parent', uitab(tabgroup, 'Title', ['Structural properties | size([x,t]) = [' num2str(length(x)) ',' num2str(length(t)) ']']));
% x_idxs = 1:min(round(modulation.lambda / diff(x(1:2))), length(x));
% t_idxs = 1:min(round(modulation.period / diff(t(1:2))), length(t));
% tile_EJ_grid = plot_structural_property(x_grid(t_idxs, x_idxs), t_grid(t_idxs, x_idxs), EJ_grid(t_idxs, x_idxs), 'Unit cell (EJ)', 'EJ [Nm^2]');
% tile_rhoA_grid = plot_structural_property(x_grid(t_idxs, x_idxs), t_grid(t_idxs, x_idxs), rhoA_grid(t_idxs, x_idxs), 'Unit cell (rhoA)', 'rhoA [Kg/m]');
% linkprop([tile_EJ_grid tile_rhoA_grid], {'View', 'XLim', 'YLim'});
% clear x_idxs t_idxs
% 
% 
% % Force and Waterfall
% axes('parent', uitab(tabgroup, 'Title', 'Force and Waterfall'));
% tile_force = plot_excitation_force(t, excitation.force);
% tile_waterfall = plot_waterfall(x, t, U, 20);


% Dispersion diagram
axes('parent', uitab(tabgroup, 'Title', ['Dispersion diagram | [Nx,Nt,Nc]=[' num2str(Nx) ',' num2str(Nt) ',' num2str(Nc) ']']));
ff_idxs = find(ff > 0 & ff < 100e3);
mu_idxs = find(mu > -5*pi & mu < 5*pi);
plot_dispersion_diagram(mu(mu_idxs), ff(ff_idxs), UU_fft(ff_idxs, mu_idxs), 'surf');
% overlay(modulation.label);
title(['Tone burst @[f,amp]=[' num2str(excitation.frequency) ',' num2str(excitation.amplitude) ']'])
ylim([0 80])