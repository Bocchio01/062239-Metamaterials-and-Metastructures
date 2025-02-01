clc
clear variables
% close all

try
    % [U, time, nodes, infos] = read_data_set('Scan_time_freezed.unv');
    [U, time, nodes, infos] = read_data_set();
catch ME
    disp([ME.message ' Interrupting.']);
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

figure('Name', ['Experimental: ' infos.filename])
colormap(jet(64))

ff_idxs = find(ff > 0 & ff < 20e3);
mu_idxs = find(mu > -5*pi & mu < 0*pi);
plot_dispersion_diagram(mu(mu_idxs), ff(ff_idxs), U_fft(ff_idxs, mu_idxs), 'surf');