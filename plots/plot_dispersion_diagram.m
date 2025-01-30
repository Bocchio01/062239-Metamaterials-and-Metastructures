function tile = plot_dispersion_diagram(mu_vector, f_vector, propagation_level, type)

arguments
    mu_vector {mustBeNumeric}
    f_vector {mustBeNumeric}
    propagation_level {mustBeNumeric} = 1
    type {mustBeText} = 'direct'
end

f_max = min(20e3, max(f_vector, [], 'all'));

switch (type)

    case 'direct'

        tile = nexttile;
        hold on
        grid on

        scatter(mu_vector / pi, f_vector * 1e-3, 5, 'or', 'filled')
        scatter(mu_vector / pi, f_vector * 1e-3, 20 * propagation_level, 'ok', 'filled')

        title('Dispersion diagram')
        xlabel('\mu/\pi')
        ylabel('f [kHz]')

        xlim([min(mu_vector) max(mu_vector)]/pi)
        ylim([0 f_max] * 1e-3)


    case 'inverse'

        tile = tiledlayout(2, 1, 'TileSpacing', 'tight');
        
        figure_alpha = nexttile;
        plot(f_vector * 1e-3, real(mu_vector), '.')
        title('Propagating part')
        xlabel('f [kHz]')
        ylabel('Re[\mu]')        
        
        figure_beta = nexttile;
        plot(f_vector * 1e-3, imag(mu_vector), '.')
        title('Attenuating part')
        xlabel('f [kHz]')
        ylabel('Im[\mu]')
        
        linkaxes([figure_alpha figure_beta], 'xy')
        xlim([0 f_max] * 1e-3)
        ylim([-pi pi]);

    otherwise
        error('Unknown dispersion type')

end

end

