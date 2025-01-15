clc
clear variables
% close all

load('parameters.mat', "beam")

%%
[U, time, nodes, infos] = read_data_set('data/nonreciprocal/Scan_time_off.unv');

dt = diff(time(1:2));
dx = diff(nodes.x(1:2));
T = diff(time([1 end]));
L = diff(nodes.x([1 end]));

t = (0:dt:T)';
x = (0:dx:L)';
U = interp2(nodes.x, time, U, x, t', "linear", 0);

Uz = zeros(size(U));
U = [
    Uz Uz Uz;
    Uz U Uz;
    Uz Uz Uz
    ];


% Frequency and wavenumber vectors
Nt = size(U, 1) - rem(size(U, 1), 2);
Nx = size(U, 2) - rem(size(U, 2), 2);

ff = (-Nt/2 : Nt/2) / (3*T);
kk = (-Nx/2 : Nx/2) * 2*pi/(3*L);
mu = kk * sum(beam.L);

% Magnitude of FFT (Power Spectrum)
fft_measured = fftshift(fft2(flipud(U)));


mLim = [-5*pi 0];
[Valmm1, pos1] = min(abs(mu-mLim(1)));
[Valmm2, pos2] = min(abs(mu-mLim(2)));
mDomain = [pos1, pos2];

fLim = [0 20e3];
[Valff1, pos1] = min(abs(ff-fLim(1)));
[Valff2, pos2] = min(abs(ff-fLim(2)));
fDomain = [pos1,pos2];

UWindow = abs(fft_measured(fDomain(1):fDomain(2), mDomain(1):mDomain(2)));
% UWindow = abs(fft_measured);
mWindow = mu(mDomain(1):mDomain(2));
fWindow = ff(fDomain(1):fDomain(2));



%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')

figure('Name', 'Experimental data')
hold on
grid on

[kk_grid, ff_grid] = meshgrid(mWindow/pi, fWindow);
contourf(kk_grid, ff_grid/1e3, UWindow/max(UWindow, [], 'all'), 0.2:0.02:1, 'EdgeColor', 'none')
% contourf(UWindow/max(UWindow, [], 'all'), 0.2:0.02:1, 'EdgeColor', 'none')
% [kk2_grid, ff2_grid] = meshgrid(mWindow/pi, fWindow);
% interp2(kk_grid, ff_grid, UWindow, )
% nexttile
% scatter(kk_grid, ff_grid/1e3)

colormap(jet(64))
colorbar
clim([0.2 1])

title('Dispersion Relation: f-k Spectrum')
xlabel('\mu/\pi [-]')
ylabel('f [kHz]')
