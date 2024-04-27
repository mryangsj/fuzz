function [S] = scattering_parallel(Z)
Q = [1, 1, 1];
Z = diag(Z);
S = 2*Q' /(Q/Z*Q') * Q/Z - eye(3);
end