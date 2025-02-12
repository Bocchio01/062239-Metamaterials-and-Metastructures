piezo  = Piezo();
Ep = piezo.E(0);


%% Instability region analysis for ideal C-, depending on C0

CN_vet = linspace(-0.5 * piezo.C_T, 1.5 * piezo.C_T, 1000);
Ep_shunted = zeros(length(CN_vet), 1);

for idx = 1:length(CN_vet)
    Ep_shunted(idx) = real(piezo.bindShunt(Shunt('RLC', 'R', 0, 'L', 0, 'C', -CN_vet(idx)), 1e10).E(0));
end


% Plot
reset(0);
set(0, 'DefaultFigureNumberTitle', 'off');
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'defaultaxesfontsize', 18);
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultTextInterpreter','latex');

figure('Name', 'Instability analysis')
tiledlayout(1, 2, 'TileSpacing', 'tight')

nexttile
hold on
grid on

plot(CN_vet / piezo.C_T, Ep_shunted / Ep, 'k', 'DisplayName', 'none', 'HandleVisibility', 'off');
xregion(-0.5, 0, 'FaceColor', 'red', 'DisplayName', 'Positive capacitance');
xregion(0, 1 - piezo.k31^2, 'FaceColor', 'green', 'DisplayName', 'Electrical instability');
yregion(-10, 0, 'FaceColor', 'blue', 'DisplayName', 'Mechanical instability');
xline(1 - piezo.k31^2, 'r--', 'HandleVisibility', 'off', 'LineWidth', 2, 'Label', '1-k_{31}^2');
yline(piezo.Y11_D / Ep, 'r--', 'HandleVisibility', 'off', 'LineWidth', 2, 'Label', 'E_p^D');

title('Stability analysis for the $C_N$ shunted piezo')
xlabel('$C_N / C_P^T$')
ylabel('$E_{SU} / E_P$')

xlim([-0.5 1.5])
ylim([-5 5])

legend('Location', 'best')

export_pdf_graphic(gcf, '/Stability_analysis_C_N')


%% Piezo Young modulus depending on frequency

shunt_labels = {'OFF', '+inf', 'C-', 'RLC-', 'RLC', 'RL//C', 'RC//L'}';

piezo  = Piezo();
Ep = piezo.E(0);

shunts = [
    % Shunt('C-')
    % Open circuit
    Shunt('+inf')

    % Short circuit
    Shunt('OFF')
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
    % Shunt('RLC', 'R', 1000, 'L', 0e-3, 'C', -30e-9);
    % Shunt('RLC', 'R', 1000, 'L', 0e-3, 'C', -7e-9);
    % Shunt('RLC', 'R', 1000, 'L', 0e-3, 'C', -5e-9);
    % 
    % % RLC-
    % % Shunt('RLC', 'R', 500, 'L', 30e-3, 'C', -30e-9);
    % % Shunt('RLC', 'R', 500, 'L', 30e-3, 'C', -7e-9);
    % % Shunt('RLC', 'R', 500, 'L', 30e-3, 'C', -5e-9);
]';



ff = logspace(0, 5, 500);
Ep_shunted = zeros(length(ff), numel(shunts));

label = cell(0);

for idx = 1:numel(shunts)

    shunt = shunts(idx);

    for ff_idx = 1:length(ff)
        Ep_shunted(ff_idx, idx)  = Piezo().bindShunt(shunt, 2*pi * ff(ff_idx)).E(0);
    end

    % label{idx} = ['L: ' num2str(shunt.L, "%.2f") 'H'];
    % label{idx} = ['RL: ' num2str(shunt.R) 'ohm' ' ' num2str(shunt.L, "%.2f") 'H'];
    label{idx} = ['RC: ' num2str(shunt.R) 'ohm' ' ' num2str(shunt.L, "%.2f") 'H' ' ' num2str(shunt.C) 'F'];

end

% Plot
reset(0);
set(0, 'DefaultFigureNumberTitle', 'off');
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'defaultaxesfontsize', 18);
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultTextInterpreter','latex');

figure('Name', 'Piezo Young modulus')
tiledlayout(2, 2, 'TileSpacing', 'tight')

tile_real = nexttile;
hold on
grid on
% xscale('log')

semilogx(ff * 1e-3, real(Ep_shunted)  / Ep);
xlim([0, 20])
ylim([-0 2])

xlabel('$f [kHz]$')
ylabel('$Re[E_{SU}] / E_p$')

legend(label, 'Location', 'best')


% tile_imag = nexttile;
% hold on
% grid on
% % xscale('log')
% 
% semilogx(ff * 1e-3, imag(Ep_shunted) / Ep)
% 
% xlim([0 20])
% title('Shunted piezo Young modulus analysis (C-, real case)')
% xlabel('$f [kHz]$')
% ylabel('$Im[E_{SU}]$')
% legend(label)
% ylim padded

% linkaxes([tile_real tile_imag], 'x')

% title('Open vs. Short shunt circuit')
% export_pdf_graphic(gcf, '/Y_SU_Open vs Short circuit')


% title('Purely inductive shunt circuit')
% export_pdf_graphic(gcf, '/Y_SU_Purely inductive shunt circuit')


% title('Resistive-Inductive shunt circuit')
% export_pdf_graphic(gcf, '/Y_SU_Resistive-Inductive shunt circuit')


% title('Resistive-(Negative) Capacitive shunt circuit')
% export_pdf_graphic(gcf, '/Y_SU_Resistive-(Negative) Capacitive shunt circuit')
