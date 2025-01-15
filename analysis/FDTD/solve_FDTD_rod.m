function [u] = solve_FDTD_rod(x_vector, t_vector, force, E_grid, rho)

u = zeros(length(t_vector), length(x_vector));
v = zeros(length(t_vector), length(x_vector));
centro = floor((length(x_vector) + 1)/2);

dx = diff(x_vector(1:2));
dt = diff(t_vector(1:2));

for q = 2:length(t_vector) - 1

    v(q, centro) = v(q, centro) + force(q);

    for m = 2:length(x_vector) - 1
        v(q+1, m) = (dt/dx)^2 / rho * E_grid(m, q) * (v(q, m+1) - 2*v(q, m) + v(q, m-1)) + 2*v(q, m) - v(q-1, m);
    end

    u(q+1, 2:length(x_vector)) = diff(v(q+1, :)) / dx;

end

end

