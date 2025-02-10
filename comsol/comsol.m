clc
clear variables
% close all

[STC, modulation] = assemble_STC('ON-OFF-OFF');


%% Parameters list

result_set = readmatrix(['comsol/results/' modulation.label '.txt']);

mu_COMSOL = result_set(:, 1) * modulation.lambda;
alpha_COMSOL = real(result_set(:, 3));
beta_COMSOL  = imag(result_set(:, 3));
polar = result_set(:, 5);

propagation_level = abs(beta_COMSOL);
propagation_level = propagation_level / norm(propagation_level, inf);



%% TMM

ff  = linspace(1, 20e3, 3000);
mu = zeros(4, length(ff));
alpha = zeros(size(mu));
beta  = zeros(size(mu));
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
    [alpha(:, ff_idx), sorting] = sort(real(mu_tmp));
    beta(:, ff_idx) = imag(mu_tmp(sorting));
    mu(:, ff_idx) = alpha(:, ff_idx) + 1i * beta(:, ff_idx);

end

clear T T1 T2 sorting ff_idx mu_tmp


%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'defaultaxesfontsize', 18);
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultTextInterpreter','latex');

figure('Name', ['/TMM_COMSOL'])
tiledlayout(2, 1)

tmm = nexttile;
hold on
grid on
grid minor

idxs = (polar >= 0.1);
scatter(alpha_COMSOL(idxs)' * 1e-3, mu_COMSOL(idxs) / pi, 20, 'ok')
scatter(alpha_COMSOL(idxs)' * 1e-3, -mu_COMSOL(idxs) / pi, 20, 'ok')

for row = 1:4
    plot(ff * 1e-3, real(mu(row, :)) / pi, '.', 'Color', [0 0.4470 0.7410]);
    plot(ff * 1e-3, imag(mu(row, :)) / pi, '.', 'Color', [0.8500 0.3250 0.0980]);
end




title('Dispersion diagram (TMM)')
ylabel('$\mu / \pi []$')
xlabel('$f [kHz]$')
ylim([-1 1]);
xlim([0 20]);

legend('Comsol solution', '', 'Propagating solutions', 'Attenuating solutions');

export_pdf_graphic(gcf, ['/TMM_COMSOL ' modulation.label ' @' num2str(modulation.omega / (2*pi) * 1e-3, 3) 'kHz'])

