function tile = plot_structural_property(x_grid, t_grid, property_grid, title_label, z_axis_label)

arguments
    x_grid {mustBeNumeric}
    t_grid {mustBeNumeric}
    property_grid {mustBeNumeric}
    title_label {mustBeText} = 'Not specified property'
    z_axis_label {mustBeText} = ''
end

tile = nexttile;
hold on
grid on
view(3)

surf(x_grid*1e3, t_grid*1e3, property_grid, 'EdgeColor', 'none')

shading interp

title(title_label)
xlabel('x [mm]')
ylabel('t [ms]')
zlabel(z_axis_label)

axis padded

end

