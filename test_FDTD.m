%% Bloch modes analysis for Discretely modulated waveguides
% Riva Emanuele - emanuele.riva@polimi.it
% Politecnico di Milano
% Department of mechanical engineering

clear
% close all
clc

%% Data

E0 = 70e9; % [Pa] Young's Modulus
rho = 2700; % [kg/m^3] density
c0 = sqrt(E0/rho);

lm = 0.072; % [m] lattice constant
am = 0.6; % [-] dimensionless modulation amplitude
nu = 0.2; % [-] dimensionless modulation speed
b = 0.02; %[m] width
vm = nu*c0;
km = 2*pi/lm;
wm = vm*km;
Tm = 2*pi/wm;

% Define sub-cells
L1 = 1/3*lm;
L2 = 2/3*lm;
L3 = 1*lm;

l(1) = L1;
l(2) = (L2 - L1);
l(3) = (L3 - L2);

Rs = 3;
% Stiffness profile
E = @(x,t) (x >= 0).*(x < L1).*E0.*(1 + am*cos( - wm*t + 0/Rs*2*pi)) + ...
           (x >= L1).*(x < L2).*E0.*(1 + am*cos( - wm*t + 1/Rs*2*pi)) + ...
           (x >= L2).*(x < L3).*E0.*(1 + am*cos( - wm*t + 2/Rs*2*pi)) + ...
           (x >= L3).*(x <= L3 + L1).*E0.*(1 + am*cos( - wm*t + 0/Rs*2*pi));

% %% Dispersion relation
% 
% P = 10; % [-] number of space harmonics
% Q = 1; % [-] number of time harmonics
% MU = 50;
% mu_min = -3*pi;
% mu_max = 3*pi;
% mu = linspace(mu_min,mu_max,MU);
% Om_min = 0;
% Om_max = 1e5;
% 
% %--- Fourier transform ---------------------------------------------------%
% 
% Npunti = 100;
% 
% x = linspace(0,lm,Npunti);
% t = linspace(0,Tm,Npunti);
% [xx,tt] = meshgrid(x,t);
% Ecell = E(xx,tt);
% 
% % Frequency domain
% ff = (-Npunti/2:Npunti/2)/Tm;
% 
% % Wavenumber domain
% Xres = 2*pi/lm;
% kkx = (-Npunti/2:Npunti/2)*Xres;
% 
% % figure(1);
% % surf(xx/lm,tt/Tm,Ecell/E0/(1+am)); shading interp; hold on;
% % colormap jet
% % cb = colorbar;
% % cb.Ticks = [1-am+0.001,1,1+am]/(1+am);
% % cb.TickLabels = {'$E_{min}$','$E_0$','$E_{max}$'};
% % cb.TickLabelInterpreter = 'LaTex';
% 
% % print(FigTag,'Fig_Cell.jpeg','-djpeg','-r600');
% 
% Ef = zeros(Npunti,Npunti);
% for ip = 1:Npunti
%     for iq = 1:Npunti
%         I = Ecell.*exp( - 1i*kkx(ip)*xx + 1i*2*pi*ff(iq)*tt);
%         em = trapz(x,I,2);
%         Ef(ip,iq) = trapz(t,em,1)/(Tm*lm);
%     end
% end
% 
% % return
% 
% %--- QEP -----------------------------------------------------------------%
% % The speed of this part can be increased a lot 
% 
% p = -P:P;
% q = -Q:Q;
% w = zeros(2*(2*Q + 1)*(2*P+1),MU);
% v = zeros((2*P+1)*(2*Q + 1),2*(2*Q + 1)*(2*P+1),MU);
% m = ones(2*(2*Q + 1)*(2*P+1),1)*mu;
% 
% for mu_idx = 1:length(mu)
%     [w(:,mu_idx), v(:,:,mu_idx)] = solve_QEP_rod(mu(mu_idx), P, Q, Ef, rho, wm, lm);
% end
% propagation_level = abs(squeeze(v(floor((2*P+1)*(2*Q+1)/2) + 1,:,:))) + 1e-10;
% 
% 
% %%
% figure
% hold on
% grid on
% 
% scatter(mu / pi, w / (2*pi), 20*propagation_level, 'ok')
% 
% title('Bloch modes')
% xlabel('\mu/\pi')
% ylabel('f [Hz]')
% 
% xlim([mu_min,mu_max]/pi);
% ylim([Om_min,Om_max])
% 
% % return

%% Numerical simulation (FDTD)

disp('FDTD')

Lf = 150*lm;
Nstep_x = 3001;
dx = Lf/(Nstep_x - 1);

Tf = Lf/c0;

xv = linspace(0,Lf,Nstep_x);
tv = linspace(0,Tf,2*Nstep_x);

% Modulation time history
EE = zeros(length(tv),length(xv));
for jj = 1:length(xv)
    EE(:,jj) = E(mod(xv(jj), lm), tv);
end


% return

%--- Input displacement --------------------------------------------------%

Om_y = 0.58;
fy = Om_y/(lm/c0);
amp = 1/(lm/c0);               % spectrum semi-width (only for cos)
nc = floor(fy*2/amp);            % number of cycles (amp = 2*fy/nc)
TIn = nc/fy;                     % duration of the force tone burst

F0 = @(tt) 1.*sin(2*pi*fy*tt).* ...
    (1-1.93*cos(2*pi*fy/nc*tt)+1.29*cos(4*pi*fy/nc*tt) ...
     -0.388*cos(6*pi*fy/nc*tt)+0.028*cos(8*pi*fy/nc*tt))/4.6360.*(tt < TIn) + ...
     +0.*((tt >= TIn));
force = F0(tv);

%--- Simulation ----------------------------------------------------------%

[UU] = solve_FDTD_rod(xv, tv, force, EE, rho);
norma = max(max(abs(UU(:,round(lm/Lf*Nstep_x):end))));

% Fourier transform - active trait

%%

L = xv(end) - xv(1);
T = Tf;

% zero padding
Uz = zeros(size(UU));
U = UU;

% Wavenumber domain
Nx = size(U,2); 
Nx = Nx - rem(Nx,2);  
Kres = 2*pi/(L); % res wavenumber [rad/m]
kk = (-Nx/2:Nx/2)*Kres;
mu = kk*lm;

% Frequency domain
Nt = size(U,1); 
Nt = Nt - rem(Nt,2); 
ff = (-Nt/2:Nt/2) * 1/(T); % res time frequency [1/s]

% FFT2D
Uf = fftshift(fft2(flipud(U/norma)))*dx*dt;


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


%% Numerical Dispersion Plot
figure

hold on
contourf(mWindow/pi, fWindow, UWindow/max(max(UWindow)),0.08:0.02:1,'EdgeColor','none');
colormap jet;
view([0,90])
clim([0.08,1])
% 
% hold on
% grid on
% 
% scatter(mu / pi, w / (2*pi), 20*propagation_level, 'ok')
% 
% title('Bloch modes')
% xlabel('\mu/\pi')
% ylabel('f [Hz]')
% xlim([mu_min,mu_max]/pi);
% ylim([Om_min,Om_max])