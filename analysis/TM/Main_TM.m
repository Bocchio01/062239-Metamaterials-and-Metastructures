clc
clear variables

[STC, modulation] = assemble_STC('ON-ON-ON');
% [STC, modulation] = assemble_STC('OFF-OFF-OFF');
STC = STC(1);


%% Dispersion relation computation via T matrix method

f  = linspace(1, 20000, 1000);
mu = zeros(4, length(f));
alpha = zeros(size(mu));
beta  = zeros(size(mu));

for f_idx = 1:length(f)

    if (strcmp(modulation.label, 'ON-ON-ON'))
        STC.Piezo = STC.Piezo.bindShunt(STC.Piezo.Shunt, 2*pi * f(f_idx));
        STC = STC.computeAveragedProps();
    end

    T1 = compute_T_beam(2*pi * f(f_idx), STC.rho{1}, STC.A{1}, STC.E{1}(0), STC.J{1}, STC.Piezo.L);
    T2 = compute_T_beam(2*pi * f(f_idx), STC.rho{2}, STC.A{2}, STC.E{2}(0), STC.J{2}, STC.Beam.L - STC.Piezo.L);

    T = T2 * T1;
    mu(:, f_idx) = sort(log(eig(T)) / 1i);
    [alpha(:, f_idx), sorting] = sort(real(mu(:, f_idx)));
    beta(:, f_idx) = imag(mu(sorting, f_idx));

end

clear T T1 T2 sorting f_idx


%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')

figure('Name', ['TM: ', modulation.label])
tiledlayout(2, 1, 'TileSpacing', 'tight')

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

plot(f*1e-3, beta, '.')
% plot(beta, f*1e-3, '.')
% plot_ScreenShot(modulation.label)

title('Attenuating part')
xlabel('f [kHz]')
ylabel('Im[\mu]')

linkaxes([figure_alpha figure_beta], 'xy')
ylim([-pi pi]);

