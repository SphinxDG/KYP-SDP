function [lambda,P,LMI,feasible] = robustSynthesisCustomSolver(A,B1,B2,C,D1,D2,Q,R,M_tilde_gen)

[n,m] = size(B1);
[d1,d2] = size(D2);
s = size(M_tilde_gen,3);

A_KYP = A';
B_KYP = [eye(n),C'];

proj1 = [zeros(n),-eye(n),zeros(n,d1);
         B1', zeros(m,n), D1'];
     
proj2 = [B2', zeros(d2,n),D2';
         zeros(d1,2*n), -eye(d1)];
     
QSR_gen = - pagemtimes(proj2',pagemtimes(M_tilde_gen,proj2));

QSR_0 = - proj1'*blkdiag(inv(Q),inv(R))*proj1;

Q_gen = QSR_0(1:n,1:n);
S_gen = QSR_0(1:n,n+1:2*n+d1);
R_gen = QSR_0(n+1:2*n+d1,n+1:2*n+d1);

Q_gen(:,:,2:s+1) = QSR_gen(1:n,1:n,:);
S_gen(:,:,2:s+1) = QSR_gen(1:n,n+1:2*n+d1,:);
R_gen(:,:,2:s+1) = QSR_gen(n+1:2*n+d1,n+1:2*n+d1,:);


N_gen = zeros(s,s,s+1);
for iii = 1:s
    N_gen(iii,iii,1+iii) = 1;
end

c = zeros(s,1);
Sigma = eye(n);

lambda = 0.01*ones(s,1);

[lambda,lambda0,~,P,Q,S,R] = customKYPLMISolver(c,Sigma,A_KYP,B_KYP,Q_gen,S_gen,R_gen,N_gen,lambda);

if lambda0 >= 0
    disp('KYP LMI infeasible')
    LMI = [];
else
    LMI = [A_KYP,B_KYP;eye(n),zeros(n,n+d1)]'*[zeros(n),P;P,zeros(n)]*[A_KYP,B_KYP;eye(n),zeros(n,n+d1)] + [Q,S;S',R];
end
feasible = lambda0<0;
end

