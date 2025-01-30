function tile = plot_fourier_coefficients(kk, ff, property_hat, title_label)

arguments
    kk {mustBeNumeric}
    ff {mustBeNumeric}
    property_hat {mustBeNumeric}
    title_label {mustBeText} = 'Not specified coefficients'
end


tile = nexttile;
hold on
grid on

imagesc(kk, ff * 1e-3, abs(property_hat));
colorbar

title(title_label)
xlabel('k [1/m]')
ylabel('f [kHz]')

axis tight

end

