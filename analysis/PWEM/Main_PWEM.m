% clc
clear variables
% close all

[STC, modulation] = assemble_STC('ON-ON-ON');


%% Structural properties evaluation

x = linspace(0, modulation.lambda, 200);
t = linspace(0, modulation.period, 20);

[x_grid, t_grid] = meshgrid(x, t);
[E_grid, J_grid, A_grid, rho_grid] = evaluate_structural_properties(x_grid, t_grid, STC);

EJ_grid = E_grid .* J_grid;
rhoA_grid = rho_grid .* A_grid;


%% Fourier decomposition
% We can't set the modulation period to 0, otherwise we can't decompose
% along the time direction the signals. We should have used a 1D FFT instead

EJ_hat = zeros(length(x), length(t));
rhoA_hat = zeros(length(x), length(t));

kk = (-length(x)/2 : length(x)/2-1) * modulation.wavenumber;
ww = (-length(t)/2 : length(t)/2-1) * 2 * pi / modulation.period;

tic
fprintf('Fourier decomposition (out of %d): %3d\n', length(kk), 0);

for kk_idx = 1:length(kk)

    for ww_idx = 1:length(ww)

        e_power = exp(-1i * (kk(kk_idx) * x_grid - ww(ww_idx) * t_grid));

        EJ_hat(kk_idx, ww_idx) = 1/modulation.period * 1/modulation.lambda * trapz(t, trapz(x, EJ_grid .* e_power, 2), 1);
        rhoA_hat(kk_idx, ww_idx) = 1/modulation.period * 1/modulation.lambda * trapz(t, trapz(x, rhoA_grid .* e_power, 2), 1);

    end

    fprintf('\b\b\b\b%3d\n', kk_idx);

end

fprintf('\b\b\b\b%.2fs\n', toc);


%% Quadratic eigenvalue problem
% The number of harmonics considered is governed by P (space) and Q (time)

mu = linspace(-5*pi, 0*pi, 100);
P = 30;
Q = 1;

alpha = zeros(2 * (2*Q+1)*(2*P+1), length(mu));
beta  = zeros((2*P+1)*(2*Q+1), 2 * (2*Q+1)*(2*P+1), length(mu));

tic
fprintf('Quadratic eigenvalue problem (out of %d): %3d\n', length(mu), 0);

for mu_idx = 1:length(mu)

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


% Structural properties
axes('parent', uitab(tabgroup, 'Title', ['Structural properties | size([x,t]) = [' num2str(length(x)) ',' num2str(length(t)) ']']));

tile_EJ_grid = plot_structural_property(x_grid, t_grid, E_grid, 'Unit cell (EJ)', 'EJ [Nm^2]');
tile_rhoA_grid = plot_structural_property(x_grid, t_grid, rhoA_grid, 'Unit cell (rhoA)', 'rhoA [Kg/m]');

linkprop([tile_EJ_grid tile_rhoA_grid], {'View', 'XLim', 'YLim'});


% Fourier coefficients
axes('parent', uitab(tabgroup, 'Title', ['Fourier coefficients | size([kk,ff]) = [' num2str(length(kk)) ',' num2str(length(ww)) ']']));

tile_EJ_hat = plot_fourier_coefficients(kk, ww / (2*pi), EJ_hat.', 'Fourier coefficients (EJ)');
tile_rhoA_hat = plot_fourier_coefficients(kk, ww / (2*pi), rhoA_hat.', 'Fourier coefficients (rhoA)');

linkaxes([tile_EJ_hat tile_rhoA_hat])
xlim([-1 1] * 2000)
ylim([-1 1] * min(max(ww / (2*pi), [], 'all'), 60) * 1e-3)


% Dispersion diagram
axes('parent', uitab(tabgroup, 'Title', ['Dispersion diagram | [P,Q]=[' num2str(P) ',' num2str(Q) ']']));
plot_dispersion_diagram(mu, alpha / (2*pi), propagation_level);
% plot_ScreenShot('ON-ON-OFF');
% ylim([6 20])
