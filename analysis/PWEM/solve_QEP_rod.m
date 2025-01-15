function [alpha, beta] = solve_QEP_rod(mu, P, Q, E_hat, rho, wm, lm)

p = -P:P;
q = -Q:Q;
km = 2*pi / lm;
k = mu / lm;

K_global = zeros((2*Q + 1)*(2*P + 1));

[SS, PP] = ndgrid(1:(2*P + 1), 1:(2*P + 1));
[RR, QQ] = ndgrid(1:(2*Q + 1), 1:(2*Q + 1));

[row_indices, col_indices] = meshgrid( ...
    (size(E_hat, 1)/2) + 1 + p(SS) - p(PP), ...
    (size(E_hat, 2)/2) + 1 + q(RR) - q(QQ) ...
    );

% Compute K_locals (4D tensor)
K_locals = E_hat(sub2ind(size(E_hat), row_indices(:), col_indices(:)));
K_locals = reshape(K_locals, 2*Q + 1, 2*Q + 1, 2*P + 1, 2*P + 1);
multiplier = (k + p * km)' .* (k + p * km);
multiplier = reshape(multiplier, 1, 1, 2*P + 1, 2*P + 1);
K_locals = K_locals .* multiplier;

% Assemble K_global
for block = 1:numel(SS)
    idxs1 = (1:2*Q + 1) + (SS(block) - 1) * (2*Q + 1);
    idxs2 = (1:2*Q + 1) + (PP(block) - 1) * (2*Q + 1);
    K_global(idxs1, idxs2) = K_global(idxs1, idxs2) + K_locals(:, :, SS(block), PP(block));
end

% Compute L* of the QEP problem
L0 = kron(eye(2*P + 1), rho*wm^2 * diag(q.^2)) - K_global;
L1 = kron(eye(2*P + 1), 2*rho*wm * diag(q));
L2 = kron(eye(2*P + 1), rho * eye(length(q)));

% Solve QEP and sort eigenvalues (and vectors)
L_norm = mean(abs([norm(L0) norm(L1) norm(L2)]));
[eig_beta, eig_alpha] = polyeig(L0/L_norm, L1/L_norm, L2/L_norm);

[alpha, sorting] = sort(real(eig_alpha));
beta = eig_beta(:, sorting);

end

