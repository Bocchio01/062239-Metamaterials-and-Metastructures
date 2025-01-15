function run_TM(modulation_model, Z_su_model, E_su_model)

load('parameters.mat'); %#ok<LOAD>
if (nargin == 0)
    [modulation, E_piezos_shunted_model, Z_su_model, E_su_model] = setup_piezos_shunted('ON-ON-ON', 'C-', 'C- (Ideal)'); %#ok<ASGLU>
else
    [modulation, E_piezos_shunted_model, Z_su_model, E_su_model] = setup_piezos_shunted(modulation_model, Z_su_model, E_su_model); %#ok<ASGLU>
end


%% Frequency domain

f = linspace(1, 20000, 1000);
w = 2*pi * f;


%% Dispersion relation computation via T matrix method

mu = zeros(4, length(f));
alpha = zeros(size(mu));
beta = zeros(size(mu));

for ii = 1:length(f)

    A_sandwich = beam.A(1) + 2*piezo.A;
    J_sandwich = beam.J(1) + 2*piezo.J;
    rho_sandwich = (beam.rho(1)*beam.A(1) + 2*piezo.rho*piezo.A) / A_sandwich;
    E_sandwich   = (beam.E(1)*beam.J(1) + 2*E_piezos_shunted_model{1}(0, w(ii))*piezo.J) / J_sandwich;

    T1 = compute_T_beam(w(ii), rho_sandwich, A_sandwich, E_sandwich, J_sandwich, beam.L(1));
    T2 = compute_T_beam(w(ii), beam.rho(2), beam.A(2), beam.E(2), beam.J(2), beam.L(2));

    T = T2 * T1;
    mu(:, ii) = sort(log(eig(T)) / 1i);
    [alpha(:, ii), sorting] = sort(real(mu(:, ii)));
    beta(:, ii) = imag(mu(sorting, ii));

end



%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')

figure('Name', 'Dispersions diagrams')
tiledlayout(1, 2)

% Dispersion relations
figure_alpha = nexttile;
hold on
grid on

plot(f*1e-3, alpha, '.')

title('Propagating part')
xlabel('f [kHz]')
ylabel('Re[\mu]')


figure_beta = nexttile;
hold on
grid on

plot(beta, f*1e-3, '.', 'DisplayName', E_su_model)

title(['Attenuating part ' modulation.model])
xlabel('f [kHz]')
ylabel('Im[\mu]')

% linkaxes([figure_alpha figure_beta])
% xlim([0 +pi])
% ylim([0 max(f*1e-3)])

legend()

linkaxes([figure_alpha figure_beta], 'xy')
% ylim([0 7])

plot_ScreenShot(modulation.model)

end