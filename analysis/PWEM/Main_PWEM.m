% clc
clear variables

% [STC, modulation] = assemble_STC('ON-ON-ON', 1);
% [STC, modulation] = assemble_STC('Sinusoidal (discrete)', 1);
[STC, modulation] = assemble_STC('Sinusoidal (continuos)', 10000);


%% Structural properties evaluation

x = linspace(0, modulation.lambda, 100);
t = linspace(0, modulation.period, 20);

[x_grid, t_grid] = meshgrid(x, t);
[E_grid, J_grid, A_grid, rho_grid] = evaluate_structural_properties(x_grid, t_grid, STC);


%% Fourier decomposition

E_hat = zeros(length(x), length(t));
J_hat = zeros(length(x), length(t));
A_hat = zeros(length(x), length(t));
rho_hat = zeros(length(x), length(t));

kk = (-length(x)/2 : length(x)/2-1) * modulation.wavenumber;
ww = (-length(t)/2 : length(t)/2-1) * modulation.omega;

tic
fprintf('Fourier decomposition (out of %d): %3d\n', length(kk), 0);

for kk_idx = 1:length(kk)

    for ww_idx = 1:length(ww)

        e_power = exp(1i * (kk(kk_idx) * x_grid - ww(ww_idx) * t_grid));

        E_hat(kk_idx, ww_idx) = 1/modulation.period * 1/modulation.lambda * trapz(t, trapz(x, E_grid .* e_power, 2), 1);
        J_hat(kk_idx, ww_idx) = 1/modulation.period * 1/modulation.lambda * trapz(t, trapz(x, J_grid .* e_power, 2), 1);
        A_hat(kk_idx, ww_idx) = 1/modulation.period * 1/modulation.lambda * trapz(t, trapz(x, A_grid .* e_power, 2), 1);
        rho_hat(kk_idx, ww_idx) = 1/modulation.period * 1/modulation.lambda * trapz(t, trapz(x, rho_grid .* e_power, 2), 1);

    end

    fprintf('\b\b\b\b%3d\n', kk_idx);

end

fprintf('Fourier decomposition done: %.2fs\n', toc)


%%

E_hat_fft = fftshift(fft2(E_grid));
J_hat_fft = fftshift(fft2(J_grid));
A_hat_fft = fftshift(fft2(A_grid));
rho_hat_fft = fftshift(fft2(rho_grid));

figure('Name', 'Fourier coefficients');

a = nexttile; imagesc(ww, kk, abs(E_hat));
b = nexttile; imagesc(ww, kk, rot90(abs(E_hat_fft)));
c = nexttile; imagesc(ww, kk, abs(J_hat));
d = nexttile; imagesc(ww, kk, rot90(abs(J_hat_fft)));

% return


%% Quadratic eigenvalue problem

mu = linspace(-3*pi, +3*pi, 100);
P = 5; % Number of harmonics in space
Q = 2; % Number of harmonics in time

alpha = zeros(2 * (2*Q+1)*(2*P+1), length(mu));
beta  = zeros((2*P+1)*(2*Q+1), 2 * (2*Q+1)*(2*P+1), length(mu));

tic
fprintf('Quadratic eigenvalue problem (out of %d): %3d\n', length(mu), 0);

for mu_idx = 1:length(mu)

    % [alpha(:, mu_idx), beta(:, :, mu_idx)] = solve_QEP_beam(mu(mu_idx), P, Q, E_hat_fft, J_hat_fft, A_hat_fft, rho_hat_fft, modulation.omega, modulation.lambda);
    % [alpha(:, mu_idx), beta(:, :, mu_idx)] = solve_QEP_beam(mu(mu_idx), P, Q, E_hat, J_hat, A_hat, rho_hat, modulation.omega, modulation.lambda);
    [alpha(:, mu_idx), beta(:, :, mu_idx)] = solve_QEP_rod_(mu(mu_idx), P, Q, E_hat, rho_hat, modulation.omega, modulation.lambda);

    fprintf('\b\b\b\b%3d\n', mu_idx);

end

fprintf('Quadratic eigenvalue problem done: %.2fs\n', toc)

propagation_level = abs(squeeze(beta(floor((2*P+1)*(2*Q+1)/2) + 1, :, :))) + 1e-10;


%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')

figure('Name', 'Generalized PWEM analysis')
tile = tiledlayout(1, 3);

plot_EJ(x_grid, t_grid, E_grid, J_grid);


nexttile(tile, 2, [1, 2])
hold on
grid on

if (modulation.omega < 2*pi*10)
    scatter(mu / pi, alpha / (2*pi) * 1e-3, 20, 'ok')
else
    scatter(mu / pi, alpha / (2*pi) * 1e-3, 20*propagation_level, 'ok')
end

title('Bloch modes')
xlabel('\mu/\pi')
ylabel('f [kHz]')

xlim([min(mu) max(mu)]/pi)
ylim([0 100])
