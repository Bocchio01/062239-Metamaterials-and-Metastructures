clc
clear variables
% close all

[STC, modulation] = assemble_STC('Sinusoidal (continuos)', 3000);


%% Numerical simulation (FDTD)

% Here E should be the mean of sandwich?
c0 = sqrt(STC(1).Beam.E / STC(1).Beam.rho);

N = 3000;
x_max = 150 * modulation.lambda;
t_max = 1 * (x_max / c0);

x = linspace(0, x_max, N);
t = linspace(0, t_max, 2*N);


% N = 3001;
% x_max = 150 * modulation.lambda;
% dx = x_max / (N-1);
% Tf = x_max/(c0)*0.6;
% 
% dt = dx/(c0*sqrt(1 + 0.6))/2;
% 
% x = linspace(0,x_max,N);
% t = linspace(0,Tf,ceil(Tf/dt));

% Modulation time history
[x_grid, t_grid] = meshgrid(x, t);
[E_grid, J_grid, A_grid, rho_grid] = evaluate_structural_properties(x_grid, t_grid, STC);


%% Excitation force
% The smaller fe or dfe, the longer the time needed

excitation.frequency = 700e+03;%8.4 * 1e3;
excitation.amplitude = 600e+03; %8.0 * 1e3;
excitation.force = tone_burst(t, excitation.frequency, excitation.amplitude);


%% Finite difference finite time

% [U] = solve_FDTD_rod(x, t, excitation.force, E_grid, rho_grid);
[U] = solve_FDTD_beam(x, t, excitation.force, E_grid, rho_grid, J_grid, A_grid);

%% Fourier transform - active trait

mu = (-length(x)/2 : length(x)/2-1) * modulation.lambda * 2*pi / x_max;
ff = (-length(t)/2 : length(t)/2-1) * 1 / t_max;

U_fft = abs(fftshift(fft2(U)));
U_fft = U_fft / max(U_fft, [], 'all');


%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')

figure('Name', ['FDTD: ' modulation.label])
tile = tiledlayout(1, 3);

tile_force = nexttile(tile, 1);
hold on
grid on

force_fft = fftshift(fft(excitation.force));
force_power_spectrum = abs(force_fft).^2;
plot(force_power_spectrum / (max(force_power_spectrum)), ff * 1e-3)

title('Force')
xlabel('A [-]')
ylabel('f [khz]')

xlim([0 1])
ylim([0 20])


% Numerical Dispersion Plot
disp_tile = nexttile(tile, 2, [1, 2]);
hold on
grid on

surf(mu / pi, ff * 1e-3, U_fft);

colormap([1 1 1; jet(64)])
colorbar
clim([0.1 1])
shading interp

title('Dispersion diagram')
xlabel('\mu / \pi [-]')
ylabel('f [khz]')

xlim([-5 5])
ylim([0 70])

linkaxes([tile_force disp_tile], 'y')