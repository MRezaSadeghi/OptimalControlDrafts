clc; close all; clear;

A = [0.5, 1; 0, 0.5];
B = [0; 1];
R = 1;
S = eye(2);
gamma = 0.95;

xdim = size(A, 1);
udim = size(B, 2);

% iteration paramere
eps_value = 1e-6;
n_iter = 10;

% initial condition
P = S;
norm_prev = norm(P);

for i = 1:n_iter
    P = S + gamma*A'*P*A  - gamma^2*A'*P*B*((R + gamma*B'*P*B)\B')*P*A;
    norm_new = norm(P);

    % iterating condition
    if abs(norm_new-norm_prev) <= eps_value
        fprintf("Converged in %d steps\n", i)
        break
    else
        norm_prev = norm_new;
    end
end

L = gamma*((R + gamma*B'*P*B)\B')*P*A;

fprintf("Final L value:\n")
disp(L)

