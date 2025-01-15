function [alpha, beta] = solve_QEP_beam(mu, P, Q, E_hat, J_hat, A_hat, rho_hat, wm, lm)

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
K_locals = -sqrt( ...
    J_hat(sub2ind(size(J_hat), row_indices(:), col_indices(:))) ./ ...
    A_hat(sub2ind(size(A_hat), row_indices(:), col_indices(:))) ...
    ) .* E_hat(sub2ind(size(E_hat), row_indices(:), col_indices(:)));
K_locals = reshape(K_locals, 2*Q + 1, 2*Q + 1, 2*P + 1, 2*P + 1);
multiplier = (k + p * km)'.^2 .* (k + p * km).^2;
multiplier = reshape(multiplier, 1, 1, 2*P + 1, 2*P + 1);
K_locals = K_locals .* multiplier;

% Assemble K_global
for block = 1:numel(SS)
    idxs1 = (1:2*Q + 1) + (SS(block) - 1) * (2*Q + 1);
    idxs2 = (1:2*Q + 1) + (PP(block) - 1) * (2*Q + 1);
    K_global(idxs1, idxs2) = K_global(idxs1, idxs2) + K_locals(:, :, SS(block), PP(block));
end

% Compute L* of the QEP problem
L0 = kron(eye(2*P + 1), rho_hat*wm^2 * diag(q.^2)) + K_global;
L1 = kron(eye(2*P + 1), 2*rho_hat*wm * diag(q));
L2 = kron(eye(2*P + 1), rho_hat * eye(length(q)));

% Solve QEP and sort eigenvalues (and vectors)
L_norm = mean(abs([norm(L0) norm(L1) norm(L2)]));
[eig_beta, eig_alpha] = polyeig(L0/L_norm, L1/L_norm, L2/L_norm);

[alpha, sorting] = sort(real(eig_alpha));
beta = eig_beta(:, sorting);

end

