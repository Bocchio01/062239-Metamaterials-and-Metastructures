function tile = plot_waterfall(x, t, U, N)

arguments
    x {mustBeNumeric}
    t {mustBeNumeric}
    U {mustBeNumeric}
    N {mustBeNumeric} = 20
end


tile = nexttile;
hold on
grid on

% U_norm = max(abs(U(:, round(size(U, 2) / 2) )), [], 'all');
U_star = U(:, round(linspace(1, size(U, 2), N)) );
x_star = x(round(linspace(1, size(U, 2), N)));
dx = 3*max(x) / N;
for ii = 1:N
    plot(t * 1e3, x_star(ii) + U_star(:, ii) / norm(U_star, inf) * dx,  'k');
end

xlabel('t [ms]')
ylabel('x [mm]')
axis tight 

end

