clc
clear variables
% close all

shunts = [
    % % Open circuit
    % Shunt('+inf')
    % 
    % % Short circuit
    % Shunt('OFF')
    % 
    % L cases
    % Shunt('RLC', 'R', 0, 'L', 0e-3, 'C', +inf);
    % Shunt('RLC', 'R', 0, 'L', 20e-3, 'C', +inf);
    % Shunt('RLC', 'R', 0, 'L', 100e-3, 'C', +inf);
    % Shunt('RLC', 'R', 0, 'L', 1000e-3, 'C', +inf);
    % 
    % % RL cases
    % Shunt('RLC', 'R', 0, 'L', 30e-3, 'C', +inf);
    % Shunt('RLC', 'R', 50, 'L', 30e-3, 'C', +inf);
    % Shunt('RLC', 'R', 200, 'L', 30e-3, 'C', +inf);
    % Shunt('RLC', 'R', 1000, 'L', 30e-3, 'C', +inf);
    % 
    % % RC-
    Shunt('RLC', 'R', 1000, 'L', 0e-3, 'C', -30e-9);
    Shunt('RLC', 'R', 1000, 'L', 0e-3, 'C', -7e-9);
    Shunt('RLC', 'R', 1000, 'L', 0e-3, 'C', -5e-9);
    % 
    % % RLC-
    % Shunt('RLC', 'R', 0, 'L', 15e-3, 'C', -100);
    % Shunt('RLC', 'R', 0, 'L', 15e-3, 'C', -30);
    % Shunt('RLC', 'R', 0, 'L', 15e-3, 'C', -5);
]';

for shunt = shunts

    [STC, modulation] = assemble_STC('ON-ON-ON', 'shunt', shunt);
    % STC = assemble_STC('ON-OFF-OFF', 'shunt', shunt);

    %% Dispersion relation computation via T matrix method

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

    figure('Name', ['/TMM_' modulation.label '_R' num2str(shunt.R) '_L' num2str(shunt.L) '_C' num2str(shunt.C)])
    tiledlayout(2, 1)

    tmm = nexttile;
    hold on
    grid on
    grid minor

    for row = 1:4
        p1 = plot(ff * 1e-3, real(mu(row, :)) / pi, '.', 'Color', [0 0.4470 0.7410]);
        p2 = plot(ff * 1e-3, imag(mu(row, :)) / pi, '.', 'Color', [0.8500 0.3250 0.0980]);
    end

    title('Dispersion diagram (TMM)')
    ylabel('$\mu / \pi []$')
    xlabel('$f [kHz]$')
    ylim([-1 1]);
    xlim([0 20]);

    legend('Propagating solutions', 'Attenuating solutions');

    % export_pdf_graphic(gcf, ['/TMM_' modulation.label '_' shunt.label '_R' num2str(shunt.R) '_L' num2str(shunt.L) '_C' num2str(shunt.C)])


    % plot_dispersion_diagram('tm (inverse) compact', mu, ff);
    % xlim([0 max(ff)] * 1e-3)

end


%% Exports
% plot_struct.export_path = 'latex/img/MATLAB';
% plot_struct.data = cell(0);
% plot_struct.data{end+1} = {figure_tag, '/title'};
% export_pdf_figure(plot_struct);