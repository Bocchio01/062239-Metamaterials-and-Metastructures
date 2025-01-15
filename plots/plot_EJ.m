function tile = plot_EJ(x_grid, t_grid, E_grid, J_grid)

tile = nexttile;
hold on
grid on
view(3)

surf(x_grid*1e3, t_grid, E_grid .* J_grid, 'EdgeColor', 'none')

title('Flexural rigidity of unit cell')
xlabel('x [mm]')
ylabel('t [s]')
zlabel('EJ [N m^2]')
axis tight

colormap(jet(64))

end

