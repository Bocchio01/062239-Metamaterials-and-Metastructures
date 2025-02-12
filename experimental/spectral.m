clc
clear variables
% close all

try
    % Load positive and negative modulations dataset
    [U_p, time_p, nodes_p, infos_p] = read_data_set();
    [U_n, time_n, nodes_n, infos_n] = read_data_set();
catch ME
    fprintf(2, [ME.message '\n']);
    return
end

%% FFT2 decomposition

N = 368 * 50;

U_fft_p = abs(fftshift(fft(U_p(:, 1), N)));
U_fft_n = abs(fftshift(fft(U_n(:, 1), N)));

Nt_p = length(U_fft_p);
Nt_n = length(U_fft_n);
ff_p = (-Nt_p/2 : Nt_p/2-1) / (Nt_p * diff(time_p(1:2)));
ff_n = (-Nt_n/2 : Nt_n/2-1) / (Nt_n * diff(time_n(1:2)));


%% Plot

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'defaultaxesfontsize', 18);
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultTextInterpreter','latex');

colormap(jet(64))
figure('Name', 'Spectral analysis')

tiledlayout(1, 2)

nexttile
hold on
grid on

idxs = find(ff_p >= 9000 & ff_p <= 12000);

norma = max([norm(U_fft_p(idxs), inf), norm(U_fft_n(idxs), inf)]);
plot(ff_p(idxs) * 1e-3, abs(U_fft_p(idxs))/norma, 'b')
plot(ff_n(idxs) * 1e-3, abs(U_fft_n(idxs))/norma, 'r')

title('Output spectra $f_e = 10.5[khz]$')
xlabel('$f [kHz]$')
ylabel('$\hat{w}_{out}$')

legend('f_m = +2[khz]', 'f_m = -2[khz]')
