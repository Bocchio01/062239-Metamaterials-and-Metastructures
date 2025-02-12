% clc
% clear variables
% close all

[STC, modulation] = assemble_STC('ON-OFF-OFF', 0);

disp('First plus, then minus')
[U1, time1, nodes1, infos1] = read_data_set();
% [U2, time2, nodes2, infos2] = read_data_set();


%% PWEM

% Structural properties evaluation

x = linspace(0, modulation.lambda, 500);
t = linspace(0, modulation.period, 100);

[x_grid, t_grid] = meshgrid(x, t);
[E_grid, J_grid, A_grid, rho_grid] = evaluate_structural_properties(x_grid, t_grid, STC);

EJ_grid = E_grid .* J_grid;
rhoA_grid = rho_grid .* A_grid;

% Fourier decomposition
fk = (-length(x)/2 : length(x)/2-1) * 1 / modulation.lambda;
ft = (-length(t)/2 : length(t)/2-1) * 1 / modulation.period;

E_hat   = 1/modulation.period * 1/modulation.lambda * fftshift(fft2(flipud(E_grid))).';
rho_hat = 1/modulation.period * 1/modulation.lambda * fftshift(fft2(flipud(rho_grid))).';
EJ_hat   = 1/modulation.period * 1/modulation.lambda * fftshift(fft2(flipud(EJ_grid))).';
rhoA_hat = 1/modulation.period * 1/modulation.lambda * fftshift(fft2(flipud(rhoA_grid))).';

% Quadratic eigenvalue problem
mu_PWEM = linspace(-5*pi, 5*pi, 200);
P = 40;
Q = 1;

alpha_PWEM = zeros(2 * (2*Q+1)*(2*P+1), length(mu_PWEM));
beta_PWEM  = zeros((2*P+1)*(2*Q+1), 2 * (2*Q+1)*(2*P+1), length(mu_PWEM));

tic
fprintf('Quadratic eigenvalue problem (out of %d): %3d\n', length(mu_PWEM), 0);

for mu_idx = 1:length(mu_PWEM)

    % [alpha(:, mu_idx), beta(:, :, mu_idx)] = solve_QEP_rod(mu(mu_idx), P, Q, E_hat, rho_hat, modulation.omega, modulation.lambda);
    [alpha_PWEM(:, mu_idx), beta_PWEM(:, :, mu_idx)] = solve_QEP_beam(mu_PWEM(mu_idx), P, Q, EJ_hat, rhoA_hat, modulation.omega, modulation.lambda);

    fprintf('\b\b\b\b%3d\n', mu_idx);

end

fprintf('\b\b\b\b%.2fs\n', toc);

propagation_level = abs(squeeze(beta_PWEM(floor((2*P+1)*(2*Q+1)/2) + 1, :, :))) + 1e-10;


%% Experimental
% 
% [U, time, nodes, infos] = read_data_set();


dt = diff(time1(1:2));
dx = mean(diff(nodes1.x));

N1 = 600;
N2 = floor(size(U1, 1) / size(U1, 2) * N1);

U_fft1 = abs(fftshift(fft2(U1, N1, N2)));
U_fft1 = fliplr(U_fft1);

% Frequency and wave-number vectors
[Nt, Nx] = size(U_fft1);
kk_EXP1 = (-Nx/2 : Nx/2-1) / (Nx * dx);
ff_EXP1 = (-Nt/2 : Nt/2-1) / (Nt * dt);

% \mu = \kappa * \lambda_m
mu_EXP1 = 2*pi*kk_EXP1 * 3*Beam().L;

% dt = diff(time2(1:2));
% dx = mean(diff(nodes2.x));
% 
% N1 = 600;
% N2 = floor(size(U2, 1) / size(U2, 2) * N1);
% 
% U_fft2 = abs(fftshift(fft2(U2, N1, N2)));
% U_fft2 = fliplr(U_fft2);
% 
% % Frequency and wave-number vectors
% [Nt, Nx] = size(U_fft2);
% kk_EXP2 = (-Nx/2 : Nx/2-1) / (Nx * dx);
% ff_EXP2 = (-Nt/2 : Nt/2-1) / (Nt * dt);
% 
% % \mu = \kappa * \lambda_m
% mu_EXP2 = 2*pi*kk_EXP2 * 3*Beam().L;


ff_idxs = find(ff_EXP1 >= 0 & ff_EXP1 <= 20e3);
mu_idxs = find(mu_EXP1 >= -5*pi & mu_EXP1 <= 0*pi);


%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'defaultaxesfontsize', 18);
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultTextInterpreter','latex');

figure_pwem = figure('Name', ['PWEM_EXP: ' modulation.label ' @' num2str(modulation.omega / (2*pi) * 1e-3, 3) 'kHz']);
colormap(jet(64))

tile = tiledlayout(2, 4);

% Structural properties
tile_EJ_grid = nexttile(tile, 1);
hold on
grid on
view(3)
surf(x_grid*1e3, t_grid, EJ_grid, 'EdgeColor', 'none')
shading interp
title('Unit cell (EJ)')
xlabel('$x [mm]$')
ylabel('$t [s]$')
zlabel('$EJ [Nm^2]$')
axis padded

tile_rhoA_grid = nexttile(tile, 5);
hold on
grid on
view(3)
surf(x_grid*1e3, t_grid, rhoA_grid, 'EdgeColor', 'none')
shading interp
title('Unit cell (rhoA)')
xlabel('$x [mm]$')
ylabel('$t [s]$')
zlabel('$rhoA [Kg/m]$')
axis padded

linkprop([tile_EJ_grid tile_rhoA_grid], {'View', 'XLim', 'YLim'});



% Dispersion diagram
dispersion_diagram = nexttile(tile, 2, [2, 3]);
hold on
grid on
contourf(mu_EXP1(mu_idxs) / pi, ff_EXP1(ff_idxs) * 1e-3, U_fft1(ff_idxs, mu_idxs), ...
    linspace(0.2, 1, 50) * max(abs(U_fft1(ff_idxs, mu_idxs)), [], 'all'), ...
    'EdgeColor', 'none', ...
    'DisplayName', 'Experimental data');
% 
% contourf(-mu_EXP2(mu_idxs) / pi, ff_EXP2(ff_idxs) * 1e-3, U_fft2(ff_idxs, mu_idxs), ...
%     linspace(0.2, 1, 50) * max(abs(U_fft2(ff_idxs, mu_idxs)), [], 'all'), ...
%     'EdgeColor', 'none', ...
%     'HandleVisibility','off');

colorbar
scatter(mu_PWEM / pi, alpha_PWEM / (2*pi)  * 1e-3, 30 * propagation_level, 'ow', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName','PWEM solution');

title('Dispersion diagram (PWEM)')
xlabel('$\mu / \pi$')
ylabel('$f [kHz]$')

xlim([min(mu_PWEM) max(mu_PWEM)]/pi)
ylim([0 20])
grid minor
legend('Experimental data', ['PWEM solution [P,Q]=[' num2str(P) ',' num2str(Q) ']'])


%%

export_pdf_graphic(gcf, ['/PWEM_EXP ' modulation.label ' @' num2str(modulation.omega / (2*pi) * 1e-3, 3) 'kHz'])
