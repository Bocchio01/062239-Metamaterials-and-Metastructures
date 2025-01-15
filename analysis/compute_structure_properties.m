function [E, J, A, rho] = compute_structure_properties(x_grid, t_grid, modulation, E_piezos_shunted_model)

load('parameters.mat'); %#ok<LOAD>

A_sandwich = beam.A(2) + 2*piezo.A;
J_sandwich = beam.J(2) + 2*piezo.J;
rho_sandwich = (beam.rho(2)*beam.A(2) + 2*piezo.rho*piezo.A) / A_sandwich;
E_sandwich = (beam.E(2)*beam.J(2) + 2*E_piezos_shunted_model*piezo.J) / J_sandwich;

x_grid = x_grid - sum(beam.L(1:6)) * floor(x_grid / sum(beam.L(1:6)));


% Models
E = distribute(modulation, E_sandwich, beam.E, x_grid, t_grid);
A = distribute(modulation, A_sandwich, beam.A, x_grid, t_grid);
J = distribute(modulation, J_sandwich, beam.J, x_grid, t_grid);
rho = distribute(modulation, rho_sandwich, beam.rho, x_grid, t_grid);

end


function properties_grid = distribute(modulation, sandwich, naked, x_grid, t_grid)

load('parameters.mat'); %#ok<LOAD>

L = cumsum(repmat([piezo.L 2], size(sandwich, 2)));

model = @(x, t) ...
    (x >= 0 & x < L(1) | x == L(6)) .* sandwich{1}(t, modulation.omega) + ...
    (x >= L(1) & x < L(2)) .* naked + ...
    (x >= L(2) & x < L(3)) .* sandwich{2}(t, modulation.omega) + ...
    (x >= L(3) & x < L(4)) .* naked + ...
    (x >= L(4) & x < L(5)) .* sandwich{3}(t, modulation.omega) + ...
    (x >= L(5) & x < L(6)) .* naked;

properties_grid = model(x_grid, t_grid);

end
