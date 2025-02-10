clc
clear variables
% close all

try
    [U1, time1, nodes1, infos1] = read_data_set();
    [U2, time2, nodes2, infos2] = read_data_set();
catch ME
    fprintf(2, [ME.message '\n']);
    return
end

%%

% N = 600;
N = 3*4891;
N = 368 * 50;

dt1 = diff(time1(1:2));
dt2 = diff(time2(1:2));

U_fft1 = abs(fftshift(fft(U1(:, 1), N)));
U_fft2 = abs(fftshift(fft(U2(:, end), N)));

Nt1 = length(U_fft1);
Nt2 = length(U_fft2);

ff_EXP1 = (-Nt1/2 : Nt1/2-1) / (Nt1 * dt1);
ff_EXP2 = (-Nt2/2 : Nt2/2-1) / (Nt2 * dt2);


% Plot

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'defaultaxesfontsize', 18);
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultTextInterpreter','latex');

colormap(jet(64))
figure()

tiledlayout(2, 2)

nexttile
hold on
grid on


idxs = find(ff_EXP1 >= 9000 & ff_EXP1 <= 12000);

norma = max([norm(U_fft1(idxs), inf), norm(U_fft2(idxs), inf)]);
plot(ff_EXP1(idxs) * 1e-3, abs(U_fft1(idxs))/norma, 'b')
plot(ff_EXP2(idxs) * 1e-3, abs(U_fft2(idxs))/norma, 'r')

% xlim([9, 12])

xlabel('$f [kHz]$')
ylabel('$\hat{w}_{out}$')

title('Output spectra $f_e = 10.5[khz]$')

legend('f_m = +2[khz]', 'f_m = -2[khz]')


% export_pdf_graphic(gcf, ['/Spectra_narrow10p5kHz_2000'])
