clc
clear variables
% close all

[STC, modulation] = assemble_STC('Sinusoidal (discrete)', 3000);


%% Preliminary analysis of wave velocities

[x_grid, t_grid] = meshgrid(linspace(0, modulation.lambda, 100), linspace(0, modulation.period, 100));
[E_grid, J_grid, A_grid, rho_grid] = evaluate_structural_properties(x_grid, t_grid, STC);

c0_longitudinal = sqrt(E_grid ./ rho_grid);
c0_transverse   = sqrt((E_grid .* J_grid) ./ (rho_grid .* A_grid));


EJ_grid = E_grid .* J_grid;
plot_structural_property(x_grid, t_grid * 1e+3, EJ_grid, '(ST) Flexural stiffness', '$EJ [Nm^2]$');

colormap(jet(64))
% export_pdf_graphic(gcf, '/ST_EJ');