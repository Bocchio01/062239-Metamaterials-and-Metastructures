function T = compute_T_beam(w, rho, A, E, J, L)

k = sqrt(w*sqrt(rho*A/(E*J)));
mu = k*L;

% T matrix coefficients
T11 = (cos(mu)+cosh(mu)) / 2;
T12 = (sin(mu)+sinh(mu)) / (2*k);
T13 = (cos(mu)-cosh(mu)) / (2*k^2*E*J);
T14 = (sin(mu)-sinh(mu)) / (2*k^3*E*J);

T21 = k/2 * (-sin(mu)+sinh(mu));
T22 = (cos(mu)+cosh(mu)) / 2;
T23 = -(sin(mu)+sinh(mu)) / (2*k*E*J);
T24 = (cos(mu)-cosh(mu)) / (2*k^2*E*J);

T31 = k^2*E*J/2 * (cos(mu)-cosh(mu));
T32 = k*E*J/2 * (sin(mu)-sinh(mu));
T33 = (cos(mu)+cosh(mu)) / 2;
T34 = (sin(mu)+sinh(mu)) / (2*k);

T41 = -k^3*E*J/2 * (sin(mu)+sinh(mu));
T42 = k^2*E*J/2 * (cos(mu)-cosh(mu));
T43 = k/2 * (-sin(mu)+sinh(mu));
T44 = 1/2 * (cos(mu)+cosh(mu));

T =  [
    T11 T12 T13 T14;
    T21 T22 T23 T24;
    T31 T32 T33 T34;
    T41 T42 T43 T44
    ];
