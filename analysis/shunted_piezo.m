
Ep = abs(Piezo().E(0));
E_su = [];

f_vet = logspace(0, 6, 500);

for f_idx = 1:length(f_vet)

    E_su(end+1)  = Piezo().bindShunt(Shunt('C-'), 2*pi*f_vet(f_idx)).E(0);

end


%%
figure

nexttile
semilogx(f_vet * 1e-3, real(E_su) / Ep)
grid on


nexttile
semilogx(f_vet * 1e-3, imag(E_su) / Ep)
grid on