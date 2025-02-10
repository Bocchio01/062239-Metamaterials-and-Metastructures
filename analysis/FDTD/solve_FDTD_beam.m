function [w] = solve_FDTD_beam(x_vector, t_vector, excitation, EJ_grid, rhoA_grid)

coordinate = floor(excitation.coordinate * length(x_vector));
offset = abs(length(x_vector) - 2*coordinate);

dx = diff(x_vector(1:2));
dt = diff(t_vector(1:2));

% Courant-Friedrichs-Lewy condition
if (max(sqrt(EJ_grid ./ rhoA_grid) * dt / dx^2, [], 'all') > 1)
    error('Stability criterion not met. Decrease dt.')
end

w = zeros(length(t_vector), offset + length(x_vector) + offset);
v = zeros(length(t_vector), offset + length(x_vector) + offset);
EJ_grid = [EJ_grid(:, -offset+end:end) EJ_grid EJ_grid(:, 1:1+offset)];
rhoA_grid = [rhoA_grid(:, -offset+end:end) rhoA_grid rhoA_grid(:, 1:1+offset)];

for t = 2:size(v, 1) - 1

    v(t, offset + coordinate) = v(t, offset + coordinate) + excitation.force(t);

    for x = 3:size(v, 2) - 2
        v(t+1, x) = -(EJ_grid(t, x) / rhoA_grid(t, x)) * (dt / dx^2)^2 * (v(t, x+2) -4*v(t, x+1) +6*v(t, x) -4*v(t, x-1) + v(t, x-2)) + 2*v(t, x) - v(t-1, x);
    end

end

x = 2:size(v, 2) - 1;
w(:, 3:end) = (v(:, x+1) -2*v(:, x) +v(:, x-1)) / dx^2;
w = w(:, offset + (1:length(x_vector)));

end

