[STC, modulation] = assemble_STC('Sinusoidal (discrete)', 'C-');


%% Numerical simulation (FDTD)

c0 = sqrt(STC(1).E{1}(0) / STC(1).rho{1});

am = 0.6;
Nx_step = 1001;

L = 100 * modulation.lambda - 5e-15;
T = 1/2 * (L / c0);
dx = L / (Nx_step - 1);
dt = dx / (2*c0*sqrt(1 + am));

x = 0:dx:L;
t = 0:dt:T;

% Modulation time history
[x_grid, t_grid] = meshgrid(x, t);
[E_grid, J_grid, A_grid, rho_grid] = evaluate_structural_properties(x_grid, t_grid, STC);
E_grid = E_grid';


%% 
%--- Input displacement --------------------------------------------------%

Om_y = 8;
fy = Om_y / (modulation.lambda/c0);
amp = 8 / (modulation.lambda/c0);               % spectrum semi-width (only for cos)
nc = floor(fy*2/amp);            % number of cycles (amp = 2*fy/nc)
TIn = nc/fy;                     % duration of the force tone burst

F0 = @(tt) 1.*sin(2*pi*fy*tt).* ...
    (1-1.93*cos(2*pi*fy/nc*tt)+1.29*cos(4*pi*fy/nc*tt) ...
     -0.388*cos(6*pi*fy/nc*tt)+0.028*cos(8*pi*fy/nc*tt))/4.6360.*(tt < TIn) + ...
     +0.*((tt >= TIn));
% F0 = @(tt) randn(length(tt), 1);
force = F0(t);


%% FDTD

% [UU] = solve_FDTD_rod(x_vector, t_vector, force, EE, beam.rho(1));
[UU] = solve_FDTD_beam(x, t, force, E_grid, rho_grid(1, :), sqrt(J_grid(1, :) ./ A_grid(1, :)));

norma = max(abs(UU(:, round(modulation.lambda/L*Nx_step):end)), [], 'all');


%% Fourier transform - active trait

% zero padding
Uz = zeros(size(UU));
U = [Uz Uz Uz;
    Uz UU Uz;
    Uz Uz Uz];

% Wavenumber domain
Nx = size(U,2); 
Nx = Nx - rem(Nx,2);  
Kres = 2*pi/(L*3); % res wavenumber [rad/m]
kk = (-Nx/2 : Nx/2) * Kres;
mu = kk * modulation.lambda;

% Frequency domain
Nt = size(U,1); 
Nt = Nt - rem(Nt,2);  
fres = 1/(T*3); % res time frequency [1/s]
ff = (-Nt/2:Nt/2)*fres;

% FFT2D
Uf = fftshift(fft2(flipud(U/norma)))*dx*dt;

mLim = [-5*pi 0*pi];
fLim = [0 20] * 1e3;

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

figure('Name', ['Modulation #' num2str(modulation.label)])
tile = tiledlayout(1, 4);

nexttile(tile, 1)

Ndisp = 20;
index = round(size(UU, 2) / Ndisp);
UVect = UU(:, round(modulation.lambda/L*Nx_step):index:end);
dL = L/Ndisp;

hold on
grid on

for ii = 1:Ndisp
    plot(t, dL*(ii-1)/modulation.lambda + dL*UVect(:,ii)/norma/2/modulation.lambda);
end

xlabel('t [ms]')
ylabel('x/\lambda_m')
ylim([0 100])


% Numerical Dispersion Plot
x_grid = mWindow/pi;
y_grid = fWindow/1e3;
[Xtmp, Ytmp] = meshgrid(x_grid, y_grid);

nexttile(tile, 2)
hold on
grid on

contourf(x_grid, y_grid, UWindow / max(UWindow, [], 'all'), 0.2:0.01:1, 'EdgeColor', 'none')

colormap(jet(64))
colorbar
clim([0.2 1])

title('Pseudo dispersion from FDTD')
xlabel('\mu / \pi [-]')
ylabel('f [khz]')



nexttile(tile);

F = abs(fft(force))/ceil(T/dt);
F = F(1:floor(ceil(T/dt)/2)+1);
F(2:end) = 2*F(2:end);
df = 1/T;
freq = (0:floor(ceil(T/dt)/2))*df;

plot(F/(max(F)), freq*(modulation.lambda / c0));

ylim([0 max(y_grid)])


% % Dispersion diagram
% nexttile
% hold on
% grid on
% 
% scatter(m(:)/pi, w(:)/2/pi*lm/c0, 10*vtmp(:), vtmp(:),'fill')
% 
% xlabel('\mu/\pi')
% ylabel('\Omega')
% xlim([mu_min mu_max]/pi);
% ylim([Om_min Om_max])

% % Force
% nexttile
% hold on
% grid on
% 
% plot(F/(max(abs(F))), freq*(modulation.lambda/c0), 'k')
% 
% xlabel('A')
% ylim([Om_min Om_max])

% nexttile
% plot_ScreenShot(modulation.model)