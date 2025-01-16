% Z_su_model = 'C-';
% [modulation, E_piezos_shunted_model, Z_su_model, E_su_model] = setup_piezos_shunted('OFF-OFF-OFF', 'Absent', 'C- (Ideal)');

% shunt = Shunt(Z_su_model);
% piezo = Piezo();
% beam  = Beam();
% STC = SpatioTemporalCell(beam, piezo);

shunt = STC.Piezo.Shunt;
piezo = STC.Piezo;
beam = STC.Beam;

%% Frequency domain

f = linspace(1, 20000, 1000);


%% Dispersion relation computation via T matrix method

mu = zeros(4, length(f));
alpha = zeros(size(mu));
beta = zeros(size(mu));

for ii = 1:length(f)

    if (~strcmp(shunt.label, 'Absent'))
        piezo = piezo.bindShunt(shunt, 2*pi*f(ii));
        STC = SpatioTemporalCell(beam, piezo);
    end

    T1 = compute_T_beam(2*pi*f(ii), STC.rho{1}, STC.A{1}, STC.E{1}(0), STC.J{1}, piezo.L);
    T2 = compute_T_beam(2*pi*f(ii), STC.rho{2}, STC.A{2}, STC.E{2}(0), STC.J{2}, beam.L - piezo.L);

    T = T2 * T1;
    mu(:, ii) = sort(log(eig(T)) / 1i);
    [alpha(:, ii), sorting] = sort(real(mu(:, ii)));
    beta(:, ii) = imag(mu(sorting, ii));

end

clear T T1 T2 sorting ii


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

plot(f*1e-3, beta, '.')

title(['Attenuating part ' shunt.label])
xlabel('f [kHz]')
ylabel('Im[\mu]')

linkaxes([figure_alpha figure_beta], 'xy')
ylim([0 7])
