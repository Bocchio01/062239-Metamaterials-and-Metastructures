clc
clear variables
% close all

[STC, modulation] = assemble_STC('ON-ON-ON');
STC = STC(1);


%% Dispersion relation computation via T matrix method

f  = linspace(1, 20e3, 1000);
mu = zeros(4, length(f));
alpha = zeros(size(mu));
beta  = zeros(size(mu));

for f_idx = 1:length(f)

    if (strcmp(modulation.label, 'ON-ON-ON'))
        STC.Piezo = STC.Piezo.bindShunt(STC.Piezo.Shunt, 2*pi * f(f_idx));
        STC = STC.computeAveragedProps();
    end

    T1 = compute_T_beam(2*pi * f(f_idx), STC.E{1}(0), STC.J{1}, STC.A{1}, STC.rho{1}, STC.Piezo.L);
    T2 = compute_T_beam(2*pi * f(f_idx), STC.E{2}(0), STC.J{2}, STC.A{2}, STC.rho{2}, STC.Beam.L - STC.Piezo.L);

    T = T2 * T1;

    mu_tmp = sort(log(eig(T)) / 1i);
    [alpha(:, f_idx), sorting] = sort(real(mu_tmp));
    beta(:, f_idx) = imag(mu_tmp(sorting));
    mu(:, f_idx) = alpha(:, f_idx) + 1i * beta(:, f_idx);

end

clear T T1 T2 sorting f_idx mu_tmp


%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')

figure('Name', ['TM: ', modulation.label])
plot_dispersion_diagram(mu, f, 0, 'inverse');
