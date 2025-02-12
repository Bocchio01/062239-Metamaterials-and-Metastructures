% clc
clear variables
% close all

[STC, modulation] = assemble_STC('Sinusoidal (continuos)', 3000);


%% Structural properties evaluation

x = linspace(0, modulation.lambda, 500);
t = linspace(0, modulation.period, 100);

[x_grid, t_grid] = meshgrid(x, t);
[E_grid, J_grid, A_grid, rho_grid] = evaluate_structural_properties(x_grid, t_grid, STC);

EJ_grid = E_grid .* J_grid;
rhoA_grid = rho_grid .* A_grid;


%% Fourier decomposition
% We can't set the modulation period to 0, otherwise we can't decompose
% along the time direction the signals. We should have used a 1D FFT instead

fk = (-length(x)/2 : length(x)/2-1) * 1 / modulation.lambda;
ft = (-length(t)/2 : length(t)/2-1) * 1 / modulation.period;

E_hat   = 1/modulation.period * 1/modulation.lambda * fftshift(fft2(flipud(E_grid))).';
rho_hat = 1/modulation.period * 1/modulation.lambda * fftshift(fft2(flipud(rho_grid))).';
EJ_hat   = 1/modulation.period * 1/modulation.lambda * fftshift(fft2(flipud(EJ_grid))).';
rhoA_hat = 1/modulation.period * 1/modulation.lambda * fftshift(fft2(flipud(rhoA_grid))).';

% EJ_hat = zeros(length(x), length(t));
% rhoA_hat = zeros(length(x), length(t));
% tic
% fprintf('Fourier decomposition (out of %d): %3d\n', length(fk), 0);
% 
% for kk_idx = 1:length(fk)
% 
%     for ff_idx = 1:length(ft)
% 
%         e_power = exp(-1i * 2*pi * (fk(kk_idx) * x_grid - ft(ff_idx) * t_grid));
% 
%         EJ_hat(kk_idx, ff_idx) = 1/modulation.period * 1/modulation.lambda * trapz(t, trapz(x, EJ_grid .* e_power, 2), 1);
%         rhoA_hat(kk_idx, ff_idx) = 1/modulation.period * 1/modulation.lambda * trapz(t, trapz(x, rhoA_grid .* e_power, 2), 1);
% 
%     end
% 
%     fprintf('\b\b\b\b%3d\n', kk_idx);
% 
% end
% 
% fprintf('\b\b\b\b%.2fs\n', toc);


%% Quadratic eigenvalue problem
% The number of harmonics considered is governed by P (space) and Q (time)

mu = linspace(-5*pi, 5*pi, 100);
P = 30;
Q = 1;

alpha = zeros(2 * (2*Q+1)*(2*P+1), length(mu));
beta  = zeros((2*P+1)*(2*Q+1), 2 * (2*Q+1)*(2*P+1), length(mu));

tic
fprintf('Quadratic eigenvalue problem (out of %d): %3d\n', length(mu), 0);

for mu_idx = 1:length(mu)

    % [alpha(:, mu_idx), beta(:, :, mu_idx)] = solve_QEP_rod(mu(mu_idx), P, Q, E_hat, rho_hat, modulation.omega, modulation.lambda);
    [alpha(:, mu_idx), beta(:, :, mu_idx)] = solve_QEP_beam(mu(mu_idx), P, Q, EJ_hat, rhoA_hat, modulation.omega, modulation.lambda);

    fprintf('\b\b\b\b%3d\n', mu_idx);

end

fprintf('\b\b\b\b%.2fs\n', toc);

propagation_level = abs(squeeze(beta(floor((2*P+1)*(2*Q+1)/2) + 1, :, :))) + 1e-10;



%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')
    
figure_pwem = figure('Name', ['PWEM: ' modulation.label ' @' num2str(modulation.omega / (2*pi) * 1e-3, 3) 'kHz']);
tabgroup = uitabgroup(figure_pwem);
colormap(jet(64))

% Create tabs
tab1 = uitab(tabgroup, 'Title', ['Structural properties | size([x,t]) = [' num2str(length(x)) ',' num2str(length(t)) ']']);
tab2 = uitab(tabgroup, 'Title', ['Fourier coefficients | size([fk,ft]) = [' num2str(length(fk)) ',' num2str(length(ft)) ']']);
tab3 = uitab(tabgroup, 'Title', ['Dispersion diagram | [P,Q]=[' num2str(P) ',' num2str(Q) ']']);


% Structural properties
axes('parent', tab1);
tile_EJ_grid = plot_structural_property(x_grid, t_grid, EJ_grid, 'Unit cell (EJ)', 'EJ [Nm^2]');
tile_rhoA_grid = plot_structural_property(x_grid, t_grid, rhoA_grid, 'Unit cell (rhoA)', 'rhoA [Kg/m]');
linkprop([tile_EJ_grid tile_rhoA_grid], {'View', 'XLim', 'YLim'});

% Fourier coefficients
axes('parent', tab2);
tile_EJ_hat = plot_fourier_coefficients(2*pi * fk, ft, EJ_hat.', 'Fourier coefficients (EJ)');
tile_rhoA_hat = plot_fourier_coefficients(2*pi * fk, ft, rhoA_hat.', 'Fourier coefficients (rhoA)');
linkaxes([tile_EJ_hat tile_rhoA_hat])
xlim([-1 1] * 2000)
ylim([-1 1] * min(max(ft, [], 'all'), 20e3) * 1e-3)

% Dispersion diagram
axes('parent', tab3);
plot_dispersion_diagram('pwem', mu, alpha / (2*pi), propagation_level);
% overlay(modulation.label, modulation.omega / (2*pi));
ylim([0 20])
grid minor

tabgroup.SelectedTab = tab3;