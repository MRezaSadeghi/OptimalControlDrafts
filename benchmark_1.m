clc; clear; close all;
% rng(5);

% dimensions
xdim = 3;
udim = 2;

% problem def
A = rand(xdim);
B = rand(xdim, udim);
S = eye(xdim);
R = eye(udim);
g = 1;

% initial value
X0 = [0.9; -0.9; 0.9];
Nk = 60; % number of iteration to converge (k)

%% Value Iteration
tic();
P = S; % initial guess

L = zeros(udim, xdim, Nk-1); % dim(L) = dimu x dimx
U = zeros(udim, 1, Nk-1);
X = zeros(xdim, 1, Nk);
Lnorm = zeros(1, Nk-1);

X(:, :, 1) = X0;

for k = 1:Nk-1
    P = S + g*A'*P*A - g^2*A'*P*B/(R + g*B'*P*B)*B'*P*A;
    L(:, :, k) = -g*(R + g*B'*P*B)\B'*P*A;
    U(:, :, k) = L(:, :, k)*X(:, :, k);
    X(:, :, k+1) = A*X(:, :, k) + B*U(:, :, k);
    Lnorm(k) = norm(L(:, :, k));
end
Lnorm_VI = Lnorm;
fprintf("VI: %2.2fms\n", toc()*1000);

%% Policy Iteration
tic();
L0 = zeros(udim, xdim); % initial guess

L = zeros(udim, xdim, Nk-1); % dim(L) = dimu x dimx
U = zeros(udim, 1, Nk-1);
X = zeros(xdim, 1, Nk);
Lnorm = zeros(1, Nk-1);

L(:, :, 1) = L0;

for k = 1:Nk-1
    Lk = L(:, :, k);
    A_coef = sqrt(g)*(A - B*Lk);
    Q_coef = S + Lk'*R*Lk;
    P = dlyap(A_coef, Q_coef);
    L(:, :, k+1) = g*(R + g*B'*P*B)\B'*P*A;
    Lnorm(k) = norm(L(:, :, k));
end

Lnorm_PI = Lnorm;
fprintf("PI: %2.2fms\n", toc()*1000);

%% Approximation
tic();
L0 = zeros(udim, xdim); % initial guess

L = zeros(udim, xdim, Nk-1); % dim(L) = dimu x dimx
U = zeros(udim, 1, Nk-1);
X = zeros(xdim, 1, Nk);
Lnorm = zeros(1, Nk-1);

X(:, :, 1) = X0;
L(:, :, 1) = L0;

r = 18;
Xs = rand(r, 3);
Phi = [Xs(:, 1).^2, Xs(:, 2).^2, Xs(:, 3).^2,...
       2*Xs(:, 1).*Xs(:, 2), 2*Xs(:, 1).*Xs(:, 3), 2*Xs(:, 2).*Xs(:, 3)];
beta = zeros(r, 1);

for k = 1:Nk-1
    Lk = L(:, :, k);
    % finding beta
    for s = 1:r
        x = Xs(s, :)';
        beta(s) = g*x'*(S + Lk'*R*Lk)*x;
        for i = 2:Nk
            x = (A + B*Lk)*x;
            beta(s) = beta(s) + g^i*x'*(S + Lk'*R*Lk)*x;    
        end
    end
    a = Phi\beta;
    P = [a(1), a(4), a(5);
         a(4), a(2), a(6);
         a(5), a(6), a(3)];
    L(:, :, k+1) = -g*(R + g*B'*P*B)\B'*P*A;
    Lnorm(k) = norm(L(:, :, k));
end

Lnorm_APP = Lnorm;
fprintf("App: %2.2fms\n", toc()*1000);

%% PLots
hold on
 
plot(1:Nk-1, Lnorm_PI, 'LineWidth', 2, "DisplayName", "PI")
plot(1:Nk-1, Lnorm_VI, 'LineWidth', 2, "DisplayName", "VI")
plot(1:Nk-1, Lnorm_APP, 'LineWidth', 2, "DisplayName", "App")
legend()