clear all
close all
clc

%Parameters list

L1 = 0.002;
L2 = 0.022;
a = (L1+L2)*3;

x1 = 0;
x2 = pi/a;

% Dispersion diagram plot

fileName = 'f_3_test2.txt'; %Test_1_no_bgs   off_off_off
A = readmatrix(fileName);
freq = A(:,3); % Frequencies [Hz]
k = A(:,1);    % Wavenumbers
freq_re = real(freq); % Real part of frequencies
polar = A(:,6); % Transversal polarization factor 

% Define the range to remove
p_min = 0.15;
p_max = 0.55;

% Create a logical index to filter out points within the specified p range

out_of_range = (polar >= p_min & polar <= p_max);

%%

%aa = scatter(k*a, freq_re/1000, 25, polar, 'filled');
aa = scatter(k*a, freq_re/1000, 25, polar, 'filled');
grid on
colorbar;
clim([0 1]); % Set color limits
xlabel('$\mu$','FontSize',16,'Interpreter','LaTex');
ylabel('f [$kHz$]','FontSize',16,'Interpreter','LaTex');
title('\textbf{Comsol}', 'Interpreter', 'latex', 'FontSize', 14);
xlim([-2*pi,2*pi])
ylim([0 20])
hold on

% Manual Bandgap definition
y_start = 0; % BG lower boundary
y_end = 0;   % BG upper boundary

% Bandgaps highlighting
bandgap = fill([0, 0, pi, pi], [y_start, y_end, y_end, y_start], [0.3010 0.7450 0.9330], 'FaceAlpha', 0.45);
hold on
bandgap.LineWidth = 1;
bandgap.EdgeColor = 'none';

y_start = 0; % BG lower boundary
y_end = 0;   % BG upper boundary

% Bandgaps highlighting
bandgap2 = fill([0, 0, pi, pi], [y_start, y_end, y_end, y_start], [0.3010 0.7450 0.9330], 'FaceAlpha', 0.45);
hold on
bandgap2.LineWidth = 1;
bandgap2.EdgeColor = 'none';

y_start = 0; % BG lower boundary
y_end = 0;   % BG upper boundary

% Bandgaps highlighting
bandgap3 = fill([0, 0, pi, pi], [y_start, y_end, y_end, y_start], [0.3010 0.7450 0.9330], 'FaceAlpha', 0.45);
hold on
bandgap3.LineWidth = 1;
bandgap3.EdgeColor = 'none';

% % Figure adjustments
% set(gcf,'Position',[100 100 800 400]) % Adjust size for two subplots in one row
% grid on
% axis square
% set(gcf, 'Color', 'White')
% set(gca,'FontSize',16)
% set(gca,'ticklabelinterpreter','LaTex')
% box on

% Adjust the colorbar and colormap
flipped_colormap = flipud(winter); 
colormap(flipped_colormap); 
c = colorbar;
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
title(c, '$p_L$', 'Interpreter', 'latex');

