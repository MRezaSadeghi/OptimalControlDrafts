clc; clear; close all;

% rng(5);

% dimensions
xdim = 2;
udim = 1;

% problem def
A = rand(xdim);
B = rand(xdim, udim);
S = eye(xdim);
R = eye(udim);
g = 1;

% initial value
X0 = [0.9; -0.9];
Nk = 60; % number of iteration to converge (k)

L = zeros(udim, xdim, Nk-1);
U = zeros(udim, 1, Nk-1);
X = zeros(xdim, 1, Nk);
Lnorm = zeros(1, Nk-1);

X(:, :, 1) = X0;
% L(:, :, 1) = L0;

%% Fitted Value Iteration (Example 11.1)

r = 6;
Us = rand(r, udim);
Xs = rand(r, xdim);
beta = zeros(r, 1);

Phi = [Xs(:, 1).^2, Xs(:, 2).^2, Us(:, 1).^2,...
       2*Xs(:, 1).*Xs(:, 2), 2*Xs(:, 1).*Us(:, 1), 2*Xs(:, 2).*Us(:, 1)];

for k = Nk:-1:1
    Gamma_k = zeros(r, 1);

    for s = 1:r
        if k == Nk
            X(:, :, Nk) = Xs(s, :)';
            xs_k = X(:, :, Nk);
            U(:, :, Nk) = Us(s, :)';
            us_k = U(:, :, Nk);
            x_plus = A*xs_k + B*us_k;
            beta(s) = x_plus'*S*x_plus;
            % beta_s at N is ready so far
        else
            X(:, :, Nk) = Xs(s, :)';
            xs_k = X(:, :, Nk);
            U(:, :, Nk) = Us(s, :)';
            us_k = U(:, :, Nk);
            x_plus = A*xs_k + B*us_k;
            beta(s) = x_plus'*(P-(rt/q)*rt')*x_plus;
        end
        Gamma_k(s, 1) = xs_k'*S*xs_k + us_k'*R*us_k+beta(s);
    end

    a = Phi\Gamma_k;
    P = [a(1), a(4); a(4), a(2)];
    rt = [a(4), a(5); a(2), a(6)];
    q = a(3);
    Lk = -q\rt';
    Lnorm(k) = norm(Lk);
    %L(:, :, k) = Lk;
end