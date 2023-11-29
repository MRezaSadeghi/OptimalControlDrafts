clc; close all; clear;

A = [0.5, 1; 0, 0.5];
B = [0; 1];
R = 1;
S = eye(2);
r = 6;
N = 10;

xdim = size(A, 1);
udim = size(B, 2);

Xs = rand(r, xdim);
Us = rand(r, udim);
beta = zeros(r, 1);

phi = [Xs(:, 1).^2, Xs(:, 2).^2, Us(:, 1).^2,...
       2*Xs(:, 1).*Xs(:, 2), 2*Xs(:, 1).*Us(:, 1), 2*Xs(:, 2).*Us(:, 1)];

L = zeros(udim, xdim, N);

for k = 1:N

    Lk = L(:, :, k);

    for s = 1:r
        xs = Xs(s, :)';
        us = Us(s, :)';

        beta(s) = xs'*S*xs + us'*R*us;
        x = A*xs + B*us;
        beta(s) = beta(s) + x'*S*x + (Lk*x)'*R*(Lk*x);

        for kk = 2:N-1
            x = (A - Lk*B)*x;
            beta(s) = beta(s) + x'*S*x + (Lk*x)'*R*(Lk*x);
        end    
    end

    a =  pinv(phi)*beta;
    P = [a(1), a(4); a(4), a(2)];
    rt = [a(5); a(6)];
    q = a(3);
    L(:, :, k+1) = q\rt';
end

L_fin = squeeze(L(:, :, N));
eig(A-B*L_fin)
L;