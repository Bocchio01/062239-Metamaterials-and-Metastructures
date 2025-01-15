function run_PWEM(modulation_model, Z_su_model, E_su_model)

load('parameters.mat'); %#ok<LOAD>
if (nargin == 0)
    [modulation, E_piezos_shunted_model, Z_su_model, E_su_model] = setup_piezos_shunted('ON-ON-ON', 'C-', 'C- (Ideal)'); %#ok<ASGLU>
else
    [modulation, E_piezos_shunted_model, Z_su_model, E_su_model] = setup_piezos_shunted(modulation_model, Z_su_model, E_su_model); %#ok<ASGLU>
end


%% Stiffness profile

x = linspace(0, modulation.lambda, 100);
t = linspace(0, modulation.period, 100);

[x_grid, t_grid] = meshgrid(x, t);
[E_grid, J_grid, A_grid, rho_grid] = compute_structure_properties(x_grid, t_grid, modulation, E_piezos_shunted_model);


%% Dispersion relation

% Fourier decomposition
E_hat = zeros(length(x), length(t));
J_hat = zeros(length(x), length(t));
A_hat = zeros(length(x), length(t));
rho_hat = zeros(length(x), length(t));

tic
for h = -length(x)/2:length(x)/2-1
    for n = -length(t)/2:length(t)/2-1

        integrand = E_grid .* exp(-1i * (h*modulation.wavenumber*x_grid - n*modulation.omega*t_grid));
        E_hat(h + length(x)/2+1, n + length(t)/2+1) = ...
            1/modulation.period * 1/modulation.lambda * trapz(t, trapz(x, integrand, 2), 1);

        integrand = J_grid .* exp(-1i * (h*modulation.wavenumber*x_grid - n*modulation.omega*t_grid));
        J_hat(h + length(x)/2+1, n + length(t)/2+1) = ...
            1/modulation.period * 1/modulation.lambda * trapz(t, trapz(x, integrand, 2), 1);

        integrand = A_grid .* exp(-1i * (h*modulation.wavenumber*x_grid - n*modulation.omega*t_grid));
        A_hat(h + length(x)/2+1, n + length(t)/2+1) = ...
            1/modulation.period * 1/modulation.lambda * trapz(t, trapz(x, integrand, 2), 1);

        integrand = rho_grid .* exp(-1i * (h*modulation.wavenumber*x_grid - n*modulation.omega*t_grid));
        rho_hat(h + length(x)/2+1, n + length(t)/2+1) = ...
            1/modulation.period * 1/modulation.lambda * trapz(t, trapz(x, integrand, 2), 1);

    end
end
toc


%% Quadratic eigenvalue problem

mu = linspace(-3*pi, 3*pi, 200);
P = 10;
Q = 1;

alpha = zeros(2 * (2*Q+1)*(2*P+1), length(mu));
beta = zeros((2*P+1)*(2*Q+1), 2 * (2*Q+1)*(2*P+1), length(mu));

for jj = 1:length(mu)
    % [alpha(:, jj), beta(:, :, jj)] = solve_QEP_rod(mu(jj), P, Q, E_hat, beam.rho(1), modulation.omega, modulation.lambda);
    [alpha(:, jj), beta(:, :, jj)] = solve_QEP_beam(mu(jj), P, Q, E_hat, J_hat, A_hat, beam.rho(1), modulation.omega, modulation.lambda);
end

propagation_level = abs(squeeze(beta(floor((2*P+1)*(2*Q+1)/2) + 1, :, :))) + 1e-10;


%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')

figure('Name', 'Generalized PWEM analysis')
tile = tiledlayout(1, 3);

plot_EJ(x_grid, t_grid, E_grid, J_grid);


nexttile(tile, [1, 2])
hold on
grid on

% propagation_level = 1;
% scatter(mu/pi, alpha / (modulation.wavenumber * sqrt(beam.E(1) / beam.rho(1))), 20*propagation_level, 'ok')
scatter(mu/pi, alpha / (modulation.wavenumber * sqrt(beam.E(1) / beam.rho(1))), 1, 'ok')

title('Bloch modes')
xlabel('\mu/\pi')
ylabel('\Omega')

xlim([-3 0])
ylim([0 3])
