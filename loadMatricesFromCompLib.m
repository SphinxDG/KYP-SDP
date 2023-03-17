function [A,B1,B2,C,D1,D2,nx,nw,nu,nz] = loadMatricesFromCompLib(name)
addpath('./Matlab_COMPlib_r1_1/COMPlib_r1_1')
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib(name);
% We consider only robust state feedback.
% To this end, we always assume C = eye(nx) and D21 = 0.


A = full(A);
B1 = full(B);
B2 = full(B);
D1 = eye(nu);
C = zeros(nu,nx);
D2 = zeros(nu,nu);
nw = nu;
nz = nu;
end

