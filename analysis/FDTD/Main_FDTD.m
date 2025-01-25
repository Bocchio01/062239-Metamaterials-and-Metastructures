clc
clear variables

[STC, modulation] = assemble_STC('Sinusoidal (continuos)', 10000);


%% Numerical simulation (FDTD)

c0 = sqrt(STC(1).Beam.E / STC(1).Beam.rho);

N = 2000;
Lf = 200 * modulation.lambda;
Tf = 1.5 * (Lf / c0);

x = linspace(0, Lf, N);
t = linspace(0, Tf, 2*N);

% Modulation time history
[x_grid, t_grid] = meshgrid(x, t);
[E_grid, J_grid, A_grid, rho_grid] = evaluate_structural_properties(x_grid, t_grid, STC);


%% Excitation force 

Om_y = 5;
fy = Om_y / (modulation.lambda/c0);
amp = 1 / (modulation.lambda/c0);               % spectrum semi-width (only for cos)
nc = floor(fy*2/amp);            % number of cycles (amp = 2*fy/nc)
TIn = nc/fy;                     % duration of the force tone burst

F0 = @(tt) 1.*sin(2*pi*fy*tt).* ...
    (1-1.93*cos(2*pi*fy/nc*tt)+1.29*cos(4*pi*fy/nc*tt) ...
     -0.388*cos(6*pi*fy/nc*tt)+0.028*cos(8*pi*fy/nc*tt))/4.6360.*(tt < TIn) + ...
     +0.*((tt >= TIn));

force = F0(t);


%% FDTD

[UU] = solve_FDTD_rod(x, t, force, E_grid, rho_grid);
% [UU] = solve_FDTD_beam(x, t, force, E_grid, rho_grid, J_grid, A_grid);
norma = max(abs(UU(:, round(modulation.lambda/Lf*N):end)), [], 'all');


%% Fourier transform - active trait

% Wavenumber domain
Nx = size(UU,2); 
Nx = Nx - rem(Nx,2);
kk = (-Nx/2 : Nx/2) * 2*pi/(Lf); % res wavenumber [rad/m]
mu = kk * modulation.lambda;

% Frequency domain
Nt = size(UU,1); 
Nt = Nt - rem(Nt,2);  
fres = 1/(Tf); % res time frequency [1/s]
ff = (-Nt/2:Nt/2)*fres;

% FFT2D
dx = diff(x(1:2));
dt = diff(t(1:2));
Uf = fftshift(fft2(flipud(UU/norma)))*dx*dt;

mLim = [-3*pi 3*pi];
fLim = [0 1e5];

[Valmm1,pos1] = min(abs(mu-mLim(1)));
[Valmm2,pos2] = min(abs(mu-mLim(2)));
mDomain = [pos1,pos2];

[Valff1,pos1] = min(abs(ff-fLim(1)));
[Valff2,pos2] = min(abs(ff-fLim(2)));
fDomain = [pos1,pos2];

UWindow = abs(Uf(fDomain(1):fDomain(2),mDomain(1):mDomain(2)));
mWindow = mu(mDomain(1):mDomain(2));
fWindow = ff(fDomain(1):fDomain(2));


%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')

figure('Name', ['FDTD: ' modulation.label])
tile = tiledlayout(1, 3);

nexttile(tile, 1)
force_fft = abs(fftshift(fft(force)));
plot(force_fft, 1:length(force_fft));

% Numerical Dispersion Plot
nexttile(tile, 2, [1, 2])
hold on
grid on

contourf(mWindow/pi, fWindow, UWindow / max(UWindow, [], 'all'), 0.08:0.01:1, 'EdgeColor', 'none')

colormap(jet(64))
colorbar
clim([0.08 1])

title('Pseudo dispersion from FDTD')
xlabel('\mu / \pi [-]')
ylabel('f [khz]')
