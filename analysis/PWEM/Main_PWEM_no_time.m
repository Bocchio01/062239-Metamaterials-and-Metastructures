% clc
clear variables

[STC, modulation] = assemble_STC('ON-ON-ON', 0);
% [STC, modulation] = assemble_STC('Sinusoidal (discrete)', 2*15000);
% [STC, modulation] = assemble_STC('Sinusoidal (continuos)', 2*15000);


%% Structural properties evaluation

x = linspace(0, modulation.lambda, 100);
t = 0;
[E_grid, J_grid, A_grid, rho_grid] = evaluate_structural_properties(x, t, STC);


%% Fourier decomposition

E_hat = zeros(length(x));
J_hat = zeros(length(x));
A_hat = zeros(length(x));
rho_hat = zeros(length(x));

kk = (-length(x)/2 : length(x)/2-1) * modulation.wavenumber;

tic
fprintf('Fourier decomposition (out of %d): %3d\n', length(kk), 0);

for kk_idx = 1:length(kk)

    e_power = exp(-1i * (kk(kk_idx) * x));
    % e_power(isnan(e_power)) = 0;

    E_hat(kk_idx, length(x)/2+1) = 1/modulation.lambda * trapz(x, E_grid .* e_power, 2);
    J_hat(kk_idx, length(x)/2+1) = 1/modulation.lambda * trapz(x, J_grid .* e_power, 2);
    A_hat(kk_idx, length(x)/2+1) = 1/modulation.lambda * trapz(x, A_grid .* e_power, 2);
    rho_hat(kk_idx, length(x)/2+1) = 1/modulation.lambda * trapz(x, rho_grid .* e_power, 2);

    fprintf('\b\b\b\b%3d\n', kk_idx);

end

fprintf('Fourier decomposition done: %.2fs\n', toc)

%%
% figure('Name', 'Fourier coefficients');
% nexttile
% imagesc(kk, ww, abs(E_hat));
% nexttile
% imagesc(kk, ww, abs(E_hat_fft));
% nexttile
% imagesc(kk, ww, abs(J_hat));
% nexttile
% imagesc(kk, ww, abs(J_hat_fft));
% E_hat = E_hat_fft * 1/modulation.period * 1/modulation.lambda;
% J_hat = J_hat_fft * 1/modulation.period * 1/modulation.lambda;
% A_hat = A_hat_fft * 1/modulation.period * 1/modulation.lambda;
% rho_hat = rho_hat_fft * 1/modulation.period * 1/modulation.lambda;


% return


%% Quadratic eigenvalue problem

mu = linspace(-3*pi, +3*pi, 50);
P = 10; % Number of harmonics in space
Q = 1; % Number of harmonics in time

alpha = zeros(2 * (2*Q+1)*(2*P+1), length(mu));
beta  = zeros((2*P+1)*(2*Q+1), 2 * (2*Q+1)*(2*P+1), length(mu));

tic
fprintf('Quadratic eigenvalue problem (out of %d): %3d\n', length(mu), 0);

for mu_idx = 1:length(mu)

    % [alpha(:, mu_idx), beta(:, :, mu_idx)] = solve_QEP_beam(mu(mu_idx), P, Q, E_hat_fft, J_hat_fft, A_hat_fft, rho_hat_fft, modulation.omega, modulation.lambda);
    [alpha(:, mu_idx), beta(:, :, mu_idx)] = solve_QEP_beam(mu(mu_idx), P, Q, E_hat, J_hat, A_hat, rho_hat, modulation.omega, modulation.lambda);
    % [alpha(:, mu_idx), beta(:, :, mu_idx)] = solve_QEP_rod_(mu(mu_idx), P, Q, E_hat, rho_hat, modulation.omega, modulation.lambda);

    fprintf('\b\b\b\b%3d\n', mu_idx);

end

fprintf('Quadratic eigenvalue problem done: %.2fs\n', toc)

propagation_level = abs(squeeze(beta(floor((2*P+1)*(2*Q+1)/2) + 1, :, :))) + 1e-10;


%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')

figure('Name', 'Generalized PWEM analysis')
% tile = tiledlayout(1, 3);
% 
% plot_EJ(x_grid, t_grid, E_grid, J_grid);


nexttile
hold on
grid on

if (modulation.omega == 0)
    scatter(mu / pi, alpha / (2*pi) * 1e-3, 20, 'ok')
else
    scatter(mu / pi, alpha / (2*pi) * 1e-3, 20*propagation_level, 'ok')
end

title('Bloch modes')
xlabel('\mu/\pi')
ylabel('f [kHz]')

xlim([min(mu) max(mu)]/pi)
ylim([0 150])
