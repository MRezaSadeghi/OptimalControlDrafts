clc; clear; close all;

A = [0.5, 1; 0, 0.5];
B = [0, 1]';

S = eye(2);
R = 1;
N = 10;

gamma = 0.95;
r = 6;
Xs = rand(r, 2);
us = rand(r, 1);
beta = zeros(r, 1);
Phi = [Xs(:, 1).^2, Xs(:, 2).^2, us.^2,...
       2*Xs(:, 1).*Xs(:, 2), 2*Xs(:, 1).*us, 2*us.*Xs(:, 2)];
X = zeros(2, N, r);
u = zeros(1, N, r);
L(1, :) = 0;
X(:, 1, :) = Xs';
u(:, 1, :) = us';

for k = 1:N
    for s = 1:r
        X = Xs(s, :)';
    
        beta(s) = X'*S*X + us'*R*us;
        for i = 1:N-1

            x_new = 
            %X(:, i+1, s) = A*X(:, i, s) + B*u(1, i, s);
            %beta(s) = beta(s) + gamma^i *X(:, i+1, s)'*S*X(:, i+1, s) + u(i)*R*u(i);
        end
        a = inv(Phi'*Phi)*Phi'*beta
        P = [a(1), a(4); a(4), a(2)];
        r = [a(4), a(5); a(2), a(6)];
        q = a(3);
    end
end
