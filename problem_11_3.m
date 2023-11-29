clc; close all; clear;

A = [0.5, 1; 0, 0.5];
B = [0; 1];
R = 1;
S = eye(2);
r = 6;
N = 10;
gamma = 0.95;

xdim = size(A, 1);
udim = size(B, 2);

% iteration paramere
eps_value = 1e-6;
n_iter = N;

% initial condition (+ samples)
rng(3);
Xs = rand(r, xdim);
Us = rand(r, udim);
L = zeros(udim, xdim, n_iter);
norm_prev = norm(L(:, :, 1));


phi = [Xs(:, 1).^2, Xs(:, 2).^2, Us(:, 1).^2,...
       2*Xs(:, 1).*Xs(:, 2), 2*Xs(:, 1).*Us(:, 1), 2*Xs(:, 2).*Us(:, 1)];


for i = 1:n_iter

    % get current Lk
    Lk = L(:, :, i);
    beta = zeros(r, 1);

    % running the model for each smaple
    x = zeros(xdim, N-1, r);

    for s = 1:r
        
        xs = Xs(s, :)';
        us = Us(s, :)';

        beta(s) = xs'*S*xs + us'*R*us;

        % find x1
        x(:, 1, s) = A*xs + B*us;

        for k = 1:N-1          
            beta(s) = beta(s) + gamma^k*(x(:, k, s)'*(S + Lk'*R*Lk)*x(:, k, s));
            x(:, k+1, s) = (A - B*Lk)*x(:, k, s);      
        end
    end

    a = phi\beta;
    % P = [a(1), a(4); a(4), a(2)];
    rt = [a(5); a(6)];
    q = a(3);
    L(:, :, i+1) = q\rt';

    norm_new = norm(L(:, :, i+1));
    % iterating condition
    if abs(norm_new-norm_prev) <= eps_value
        fprintf("Converged in %d steps\n", i)
        break
    else
        norm_prev = norm_new;
    end

end

figure(1)
for i = 1:r
    x1 = x(1, :, i);
    x2 = x(2, :, i);
    plot(x1)
    plot(x2, '--')
    hold on
end

L_fin = squeeze(L(:, :, i-1));
fprintf("Final L value:\n")
disp(L_fin)