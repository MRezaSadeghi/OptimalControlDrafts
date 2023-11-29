clc; clear; close all;

A = [0.5, 1; 0, 0.5];
B = [0; 1];
R = 1;
S = eye(2);

gamma = 0.95;
e0 = 1e-6;

n = 10;
L = zeros(n, 2);


L(1, :) = 0;
e = 1;
for i = 1:10
    Lmat = L(i, :);
    %Lmat'*R*Lmat
    %S + Lmat'*R*Lmat
    P = dlyap(sqrt(gamma)*(A - B*Lmat)', S + Lmat'*R*Lmat);
    L(i+1, :) = gamma*(R + gamma*B'*P*B)\(B'*P*A);
    delta = norm(L(i+1, :) - L(i, :));
    if delta < e0
        break;
    end
end


