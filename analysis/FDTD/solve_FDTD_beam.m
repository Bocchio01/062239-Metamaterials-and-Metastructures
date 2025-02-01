function [w] = solve_FDTD_beam(x_vector, t_vector, excitation, E_grid, J_grid, A_grid, rho_grid)

coordinate = floor(excitation.coordinate * length(x_vector));
offset = abs(length(x_vector) - 2*coordinate);

dx = diff(x_vector(1:2));
dt = diff(t_vector(1:2));

w = zeros(length(t_vector), offset + length(x_vector) + offset);
v = zeros(length(t_vector), offset + length(x_vector) + offset);
E_grid = [E_grid(:, -offset+end:end) E_grid E_grid(:, 1:1+offset)];
J_grid = [J_grid(:, -offset+end:end) J_grid J_grid(:, 1:1+offset)];
A_grid = [A_grid(:, -offset+end:end) A_grid A_grid(:, 1:1+offset)];
rho_grid = [rho_grid(:, -offset+end:end) rho_grid rho_grid(:, 1:1+offset)];

for t = 2:size(v, 1) - 1

    v(t, offset + coordinate) = v(t, offset + coordinate) + excitation.force(t);

    for x = 3:size(v, 2) - 2
        v(t+1, x) = (dt^2 / -rho_grid(t, x)) * (J_grid(t, x) / A_grid(t, x) * E_grid(t, x) * ( v(t, x+2) - 4*v(t, x+1) + 6*v(t, x) -4*v(t, x-1) + v(t, x-2) ) / dx^4) + 2*v(t, x) - v(t-1, x);
    end

    w(t+1, 3:end) = diff(diff(v(t+1, :))) / dx^2;

    if (any(isnan(w), 'all'))
        error(['Simulation diverged @t=' num2str(t)])
    end

end

end

