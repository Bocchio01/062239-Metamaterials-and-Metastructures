clc
clear variables
% close all

try
    % [U, time, nodes, infos] = read_data_set('Scan_time_freezed.unv');
    [U, time, nodes, infos] = read_data_set();
catch ME
    fprintf(2, [ME.message '\n']);
    return
end


dt = diff(time(1:2));
dx = mean(diff(nodes.x));

N1 = 600;
N2 = floor(size(U, 1) / size(U, 2) * N1);

U_fft = abs(fftshift(fft2(U, N1, N2)));
U_fft = fliplr(U_fft);

% Frequency and wave-number vectors
[Nt, Nx] = size(U_fft);
kk = (-Nx/2 : Nx/2-1) / (Nx * dx);
ff = (-Nt/2 : Nt/2-1) / (Nt * dt);

% \mu = \kappa * \lambda_m
mu = 2*pi*kk * 3*Beam().L;


%% Plot

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off')
set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'defaultaxesfontsize', 18);
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultTextInterpreter','latex');

figure('Name', ['Experimental: ' infos.filename])

colormap(jet(64))

tiledlayout(1, 2)

ff_idxs = find(ff >= 0 & ff <= 20e3);
mu_idxs = find(mu >= -5*pi & mu <= 0*pi);
plot_dispersion_diagram('experimental', -mu(mu_idxs), ff(ff_idxs), U_fft(ff_idxs, mu_idxs));

title('Experimetal fft2()')
xlabel('$\mu / \pi$')
ylabel('$f [kHz]$')

ylim([0 20])

idxs = find(time > 0e-3);
plot_waterfall(flip(nodes.x) * 1e+3, time(idxs), U(idxs, :), 15);

title('Waterfall plot')
