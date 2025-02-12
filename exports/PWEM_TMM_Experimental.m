clc
clear variables
% close all

[STC, modulation] = assemble_STC('ON-OFF-OFF', 0);

[U, time, nodes, infos] = read_data_set();

%% TMM

ff  = linspace(1, 20e3, 1000);
mu_TMM = zeros(4, length(ff));
alpha_TMM = zeros(size(mu_TMM));
beta_TMM  = zeros(size(mu_TMM));
TT = zeros(4, 4, 2*numel(STC));

for ff_idx = 1:length(ff)

    for STC_idx = 1:numel(STC)

        % Current STC
        cSTC = STC(STC_idx);

        cSTC.Piezo = cSTC.Piezo.bindShunt(cSTC.Piezo.Shunt, 2*pi * ff(ff_idx));
        cSTC = cSTC.computeAveragedProps();

        TT(:, :, 2*(STC_idx-1)+1) = compute_T_beam(2*pi * ff(ff_idx), cSTC.E{1}(0), cSTC.J{1}, cSTC.A{1}, cSTC.rho{1}, cSTC.Piezo.L);
        TT(:, :, 2*(STC_idx-1)+2) = compute_T_beam(2*pi * ff(ff_idx), cSTC.E{2}(0), cSTC.J{2}, cSTC.A{2}, cSTC.rho{2}, cSTC.Beam.L - cSTC.Piezo.L);

    end

    T = TT(:, :, end);
    for k = size(TT, 3)-1:-1:1
        T = T * TT(:, :, k);
    end

    mu_tmp = sort(log(eig(T)) / 1i);
    [alpha_TMM(:, ff_idx), sorting] = sort(real(mu_tmp));
    beta_TMM(:, ff_idx) = imag(mu_tmp(sorting));
    mu_TMM(:, ff_idx) = alpha_TMM(:, ff_idx) + 1i * beta_TMM(:, ff_idx);

end


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
mu_PWEM = linspace(-5*pi, -0*pi, 100);
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


dt = diff(time(1:2));
dx = mean(diff(nodes.x));

N1 = 600;
N2 = floor(size(U, 1) / size(U, 2) * N1);

U_fft = abs(fftshift(fft2(U, N1, N2)));
U_fft = fliplr(U_fft);

% Frequency and wave-number vectors
[Nt, Nx] = size(U_fft);
kk_EXP = (-Nx/2 : Nx/2-1) / (Nx * dx);
ff_EXP = (-Nt/2 : Nt/2-1) / (Nt * dt);

% \mu = \kappa * \lambda_m
mu_EXP = 2*pi*kk_EXP * 3*Beam().L;


%% Plot

figure('Name', ['Experimental: ' infos.filename])
colormap(jet(64))

ff_idxs = find(ff_EXP >= 0 & ff_EXP <= 20e3);
mu_idxs = find(mu_EXP >= -5*pi & mu_EXP <= 0*pi);
plot_dispersion_diagram('experimental', mu_EXP(mu_idxs), ff_EXP(ff_idxs), U_fft(ff_idxs, mu_idxs));

idxs = find(time > 3e-3);
plot_waterfall(flip(nodes.x), time(idxs), U(idxs, :), 16);


%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'defaultaxesfontsize', 18);
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultTextInterpreter','latex');

figure_pwem = figure('Name', ['PWEM_TMM_EXP: ' modulation.label ' @' num2str(modulation.omega / (2*pi) * 1e-3, 3) 'kHz']);
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
dispersion_diagram = nexttile(tile, 2, [2, 2]);
hold on
grid on
contourf(mu_EXP(mu_idxs) / pi, ff_EXP(ff_idxs) * 1e-3, U_fft(ff_idxs, mu_idxs) / max(abs(U_fft(ff_idxs, mu_idxs)), [], 'all'), ...
    linspace(0.2, 1, 50), ...
    'EdgeColor', 'none', ...
    'DisplayName', 'Experimental data');

colorbar
scatter(mu_PWEM / pi, alpha_PWEM / (2*pi)  * 1e-3, 30 * propagation_level, 'ow', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName','PWEM solution');

title('Dispersion diagram (PWEM)')
xlabel('$\mu / \pi$')
ylabel('$f [kHz]$')

xlim([min(mu_PWEM) max(mu_PWEM)]/pi)
ylim([0 20])
grid minor
legend('Experimental data', ['PWEM solution [P,Q]=[' num2str(P) ',' num2str(Q) ']'])


result_set = readmatrix(['comsol/results/' modulation.label '.txt']);

mu_COMSOL = result_set(:, 1) * modulation.lambda;
alpha_COMSOL = real(result_set(:, 3));
beta_COMSOL  = imag(result_set(:, 3));
polar = result_set(:, 5);

propagation_level_COMSOL = abs(beta_COMSOL);
propagation_level_COMSOL = propagation_level_COMSOL / norm(propagation_level_COMSOL, inf);

% TMM
tmm = nexttile(tile, 4, [2, 1]);
hold on
grid on
grid minor

idxs = (polar >= 0.1);
scatter(mu_COMSOL(idxs) / pi, alpha_COMSOL(idxs) * 1e-3, 20, 'ok')
scatter(-mu_COMSOL(idxs) / pi, alpha_COMSOL(idxs) * 1e-3, 20, 'ok')

for row = 1:4
    p1 = plot(real(mu_TMM(row, :)) / pi, ff * 1e-3, '.', 'Color', [0 0.4470 0.7410]);
    p2 = plot(imag(mu_TMM(row, :)) / pi, ff * 1e-3, '.', 'Color', [0.8500 0.3250 0.0980]);
end

title('Dispersion diagram (TMM)')
xlabel('$\mu / \pi []$')
ylabel('$f [kHz]$')
xlim([-1 1]);
ylim([0 20]);

legend('Comsol solution', '', 'Propagating solutions', 'Attenuating solutions');


%%

export_pdf_graphic(gcf, ['/PWEM_TMM_EXP ' modulation.label ' @' num2str(modulation.omega / (2*pi) * 1e-3, 3) 'kHz'])
