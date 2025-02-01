function tile = plot_excitation_force(t, force)

arguments
    t {mustBeNumeric}
    force {mustBeNumeric}
end

tile = nexttile;
hold on
grid on

plot(t, force)

title('Excitation force')
xlabel('t [s]')
ylabel('f [N]')
axis tight


tile = nexttile;
hold on
grid on

force_fft = fftshift(fft(force));
force_power_spectrum = abs(force_fft).^2;
ff = (-length(t)/2 : length(t)/2-1) * 1 / max(t);
plot(ff * 1e-3, force_power_spectrum, 'o')

title('Excitation force (fft)')
xlabel('f [khz]')
ylabel('A [-]')

xlim([0 100])
ylim([0 max(force_power_spectrum)])

end

