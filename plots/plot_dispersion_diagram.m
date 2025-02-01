function tile = plot_dispersion_diagram(mu_vector, f_vector, propagation_level, type, N_levels)

arguments
    mu_vector {mustBeNumeric}
    f_vector {mustBeNumeric}
    propagation_level {mustBeNumeric} = 1
    type {mustBeText} = 'dotted'
    N_levels {mustBeNumeric} = 100
end

f_max = min(20e3, max(f_vector, [], 'all'));

switch (type)

    case 'dotted'

        if (all(isreal(mu_vector)))

            % PWEM case
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

        else

            % TM case
            tile = tiledlayout(1, 2, 'TileSpacing', 'tight');
            
            figure_alpha = nexttile;
            hold on
            grid on
            plot(real(mu_vector), f_vector * 1e-3, '.')
            title('Propagating part')
            ylabel('f [kHz]')
            xlabel('Re[\mu]')        
            
            figure_beta = nexttile;
            hold on
            grid on
            plot(imag(mu_vector), f_vector * 1e-3, '.')
            title('Attenuating part')
            ylabel('f [kHz]')
            xlabel('Im[\mu]')
            
            linkaxes([figure_alpha figure_beta], 'xy')
            ylim([0 f_max] * 1e-3)
            xlim([-pi pi]);
        
        end


    case 'surf'

        tile = nexttile;
        hold on
        grid on
        
        contourf(mu_vector / pi, f_vector * 1e-3, propagation_level, ...
            linspace(0.2, 1, N_levels) * max(abs(propagation_level), [], 'all'), ...
            'EdgeColor', 'none');

        colorbar

        xlim([min(mu_vector) max(mu_vector)]/pi)
        ylim([0 f_max] * 1e-3)


    otherwise
        error('Unknown dispersion type')

end

end

