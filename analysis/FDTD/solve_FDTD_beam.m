function [w] = solve_FDTD_beam(x_vector, t_vector, force, E_grid, rho, Rg)

w = zeros(length(t_vector), length(x_vector));
v = zeros(length(t_vector), length(x_vector));
centro = floor((length(x_vector) + 1)/2);

dx = diff(x_vector(1:2));
dt = diff(t_vector(1:2));

for q = 2:length(t_vector) - 1

    v(q, centro) = v(q, centro) + force(q);
    
    for m = 3:length(x_vector) - 2
         v(q+1, m) = (dt^2 / -rho(m)) * (Rg(m)^2 * E_grid(m, q) * ( v(q, m+2) - 4*v(q, m+1) + 6*v(q, m) -4*v(q, m-1) + v(q, m-2) ) / dx^4) + 2*v(q, m) - v(q-1, m);
    end
    
    w(q+1, 2:length(x_vector)-1) = diff(v(q+1, :), 2) / dx^2;

end

end

