% Erivelton Gualter
%
% HW3


% Problem 4

A = [0 1; -12 -8];

P = lyap(A, eye(2))

syms x1 x2

Vx = [x1 x2]*P*[x1; x2]

pretty(simplify(Vx))

Q = -(A'*P+P*A)

eig(Q)

