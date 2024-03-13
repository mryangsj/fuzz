function [S] = scattering_series(Z)
B = [1, 1, 1];
Z = diag(Z);
S = eye(3) - 2*Z*B' /(B*Z*B') * B;
end