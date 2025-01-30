function [alpha, beta] = solve_QEP_beam(mu, P, Q, EJ_hat, rhoA_hat, wm, lm)

p = -P:P;
q = -Q:Q;
km = 2*pi / lm;
k = mu / lm;

K_global = zeros((2*Q + 1)*(2*P + 1));
M0_global = zeros((2*Q + 1)*(2*P + 1));
M1_global = zeros((2*Q + 1)*(2*P + 1));
M2_global = zeros((2*Q + 1)*(2*P + 1));

[SS, PP] = ndgrid(1:(2*P + 1), 1:(2*P + 1));
[RR, QQ] = ndgrid(1:(2*Q + 1), 1:(2*Q + 1));

[row_indices, col_indices] = meshgrid( ...
    (size(EJ_hat, 1)/2) + 1 + p(SS) - p(PP), ...
    (size(EJ_hat, 2)/2) + 1 + q(RR) - q(QQ) ...
    );

% Compute K_locals (4D tensor)
K_locals = EJ_hat(sub2ind(size(EJ_hat), row_indices(:), col_indices(:)));
[M0_locals, M1_locals, M2_locals] = deal(rhoA_hat(sub2ind(size(rhoA_hat), row_indices(:), col_indices(:))));

K_locals = reshape(K_locals, 2*Q + 1, 2*Q + 1, 2*P + 1, 2*P + 1);
M0_locals = reshape(M0_locals, 2*Q + 1, 2*Q + 1, 2*P + 1, 2*P + 1) .* (meshgrid(q) .* meshgrid(q)') * wm^2;
M1_locals = reshape(M1_locals, 2*Q + 1, 2*Q + 1, 2*P + 1, 2*P + 1) .* (meshgrid(q) + meshgrid(q)') * wm;
M2_locals = reshape(M2_locals, 2*Q + 1, 2*Q + 1, 2*P + 1, 2*P + 1);

K_locals_multiplier = (k + p * km)'.^2 .* (k + p * km).^2;
K_locals_multiplier = reshape(K_locals_multiplier, 1, 1, 2*P + 1, 2*P + 1);

K_locals = K_locals .* K_locals_multiplier;

% Assembly global matrices
for block = 1:numel(SS)
    
    idxs1 = (1:2*Q + 1) + (SS(block) - 1) * (2*Q + 1);
    idxs2 = (1:2*Q + 1) + (PP(block) - 1) * (2*Q + 1);

    K_global(idxs1, idxs2) = K_global(idxs1, idxs2) + K_locals(:, :, SS(block), PP(block));
    M0_global(idxs1, idxs2) = M0_global(idxs1, idxs2) + M0_locals(:, :, SS(block), PP(block));
    M1_global(idxs1, idxs2) = M1_global(idxs1, idxs2) + M1_locals(:, :, SS(block), PP(block));
    M2_global(idxs1, idxs2) = M2_global(idxs1, idxs2) + M2_locals(:, :, SS(block), PP(block));

end

% Compute L* of the QEP problem
L0 = -K_global + M0_global;
L1 = M1_global;
L2 = M2_global;

% Solve QEP and sort eigenvalues (and vectors)
L_norm = mean(abs([norm(L0) norm(L1) norm(L2)]));
[eig_beta, eig_alpha] = polyeig(L0/L_norm, L1/L_norm, L2/L_norm);

[alpha, sorting] = sort(real(eig_alpha));
beta = eig_beta(:, sorting);

end

