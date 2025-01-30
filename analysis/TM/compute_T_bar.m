function T = compute_T_bar(w, E, A, rho, L)

c  = sqrt(E / rho);
k  = w / c;
mu = k * L;

% T matrix coefficients
T11 = cos(mu);
T12 = 1 / (E*A*k) * sin(mu);

T21 = -E*A*k * sin(mu);
T22 = cos(mu);

T = [
    T11 T12;
    T21 T22
    ];
