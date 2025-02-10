function tile = plot_dispersion_diagram(type, mu_vector, f_vector, data, N_levels)

arguments
    type {mustBeText}
    mu_vector {mustBeNumeric}
    f_vector {mustBeNumeric}
    data {mustBeNumeric} = 1
    N_levels {mustBeNumeric} = 100
end

f_max = min(20e3, max(f_vector, [], 'all'));

switch (lower(type))

    case 'pwem'

        tile = nexttile;
        hold on
        grid on

        scatter(mu_vector / pi, f_vector * 1e-3, 5, 'or', 'filled')
        scatter(mu_vector / pi, f_vector * 1e-3, 20 * data, 'ok', 'filled')

        title('Dispersion diagram')
        xlabel('\mu/\pi')
        ylabel('f [kHz]')

        xlim([min(mu_vector) max(mu_vector)]/pi)
        ylim([0 f_max] * 1e-3)

    case 'tm (inverse)'

        tile = tiledlayout(2, 1, 'TileSpacing', 'tight');

        figure_alpha = nexttile;
        hold on
        grid on
        grid minor
        plot(f_vector * 1e-3, real(mu_vector), '.')
        title('Propagating part')
        xlabel('f [kHz]')
        ylabel('Re[\mu]')

        figure_beta = nexttile;
        hold on
        grid on
        grid minor
        plot(f_vector * 1e-3, imag(mu_vector), '.')
        title('Attenuating part')
        xlabel('f [kHz]')
        ylabel('Im[\mu]')

        linkaxes([figure_alpha figure_beta], 'xy')
        xlim([0 f_max] * 1e-3)
        ylim([-pi pi]);

    case 'tm (direct)'

        tile = tiledlayout(1, 2, 'TileSpacing', 'tight');

        figure_alpha = nexttile;
        hold on
        grid on
        grid minor
        plot(real(mu_vector), f_vector * 1e-3, '.')
        title('Propagating part')
        xlabel('Re[\mu]')
        ylabel('f [kHz]')

        figure_beta = nexttile;
        hold on
        grid on
        grid minor
        plot(imag(mu_vector), f_vector * 1e-3, '.')
        title('Attenuating part')
        xlabel('Im[\mu]')
        ylabel('f [kHz]')

        linkaxes([figure_alpha figure_beta], 'xy')
        xlim([-pi pi]);
        ylim([0 f_max] * 1e-3)
        
    case 'tm (inverse) compact'

        tile = nexttile;
        hold on
        grid on
        grid minor
        plot(f_vector * 1e-3, real(mu_vector), 'k.', 'DisplayName', 'Propagating solutions')
        plot(f_vector * 1e-3, imag(mu_vector), 'r.', 'DisplayName', 'Attenuating solutions')
        title('??')
        xlabel('f [kHz]')
        ylabel('Im[\mu]')

        xlim([0 f_max] * 1e-3)
        ylim([0 pi]);


    case {'fdtd', 'experimental'}

        tile = nexttile;
        hold on
        grid on

        contourf(mu_vector / pi, f_vector * 1e-3, data, ...
            linspace(0.2, 1, N_levels) * max(abs(data), [], 'all'), ...
            'EdgeColor', 'none');

        colorbar

        title('Dispersion diagram')
        xlim([min(mu_vector) max(mu_vector)]/pi)
        ylim([0 f_max] * 1e-3)

    otherwise
        error('Unknown dispersion type')

end

end

