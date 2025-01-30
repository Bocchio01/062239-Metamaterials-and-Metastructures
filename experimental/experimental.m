clc
clear variables
% close all

[U, time, nodes, infos] = read_data_set('Scan_time_freezed.unv');
% [U, time, nodes, infos] = read_data_set();

dt = diff(time(1:2));
dx = mean(diff(nodes.x));

N1 = 600;
N2 = floor(size(U, 1) / size(U, 2) * N1);

U_fft = abs(fftshift(fft2(U, N1, N2)));
U_fft = fliplr(U_fft);
U_fft = U_fft / max(U_fft, [], 'all');

% Frequency and wave-number vectors
[Nt, Nx] = size(U_fft);
kk = (-Nx/2 : Nx/2-1) / (Nx * dx);
ff = (-Nt/2 : Nt/2-1) / (Nt * dt);

% \mu = \kappa * \lambda_m
mu = 2*pi*kk * 3*Beam().L;


%% Plot

figure('Name', ['Experimental: ' infos.filename])

nexttile
hold on
grid on

surf(mu / pi, ff * 1e-3, U_fft);
alpha 0.5
% plot_ScreenShot('ON-ON-ON')
% contour(mu/pi, ff * 1e-3, U_power_spectra, 0.02:0.1:1);
% contourf(mu/pi, ff * 1e-3, U_power_spectra, 0.02:0.01:1, 'EdgeColor', 'none');
% pcolor(mu/pi, ff * 1e-3, U_power_spectra);
% imagesc(mu/pi, ff * 1e-3, U_power_spectra);

colormap([1 1 1; jet(64)])
colorbar
clim([0.2 1])
shading interp

xlim([-5 0])
ylim([0 20])

% xlim([-5 -3.5])
% ylim([6 13])