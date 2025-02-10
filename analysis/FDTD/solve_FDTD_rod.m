function [u] = solve_FDTD_rod(x_vector, t_vector, excitation, E_grid, rho_grid)

coordinate = floor(excitation.coordinate * length(x_vector));
offset = abs(length(x_vector) - 2*coordinate);

dx = diff(x_vector(1:2));
dt = diff(t_vector(1:2));

% Courant-Friedrichs-Lewy condition
if (max(sqrt(E_grid ./ rho_grid) * dt/dx, [], 'all') > 1)
    error('Stability criterion not met. Decrease dt.')
end

u = zeros(length(t_vector), offset + length(x_vector) + offset);
v = zeros(length(t_vector), offset + length(x_vector) + offset);
E_grid = [E_grid(:, -offset+end:end) E_grid E_grid(:, 1:1+offset)];
rho_grid = [rho_grid(:, -offset+end:end) rho_grid rho_grid(:, 1:1+offset)];

for t = 2:size(v, 1) - 1

    v(t, offset + coordinate) = v(t, offset + coordinate) + excitation.force(t);

    for x = 2:size(v, 2) - 1
        v(t+1, x) = (dt/dx)^2 / rho_grid(t, x) * E_grid(t, x) * (v(t, x+1) - 2*v(t, x) + v(t, x-1)) + 2*v(t, x) - v(t-1, x);
    end

end

x = 2:size(v, 2) - 1;
u(:, 3:end) = (v(:, x+1) -v(:, x-1)) / dx;
u = u(:, offset + (1:length(x_vector)));

end

