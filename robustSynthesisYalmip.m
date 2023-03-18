function [lambda,P,LMI,feasible] = robustSynthesisYalmip(A,B1,B2,C,D1,D2,Q,R,M_tilde_gen,solver)

[n,m] = size(B1);
[d1,d2] = size(D2);
s = size(M_tilde_gen,3);
proj0 = [-eye(n), zeros(n,n+d1);
         A', eye(n), C';
         zeros(n),-eye(n),zeros(n,d1);
         B1', zeros(m,n), D1'
         B2', zeros(d2,n),D2';
         zeros(d1,2*n), -eye(d1)];
     
P = sdpvar(n,n);
lambda = sdpvar(s,1);

M_tilde = zeros(size(M_tilde_gen,1),size(M_tilde_gen,1),'like',sdpvar);
for iii = 1:s
    M_tilde = M_tilde + lambda(iii)*M_tilde_gen(:,:,iii);
end

LMI = proj0'*blkdiag([zeros(n),P;P,zeros(n)],inv(Q),inv(R),M_tilde)*proj0;

ops = sdpsettings('solver',solver,'verbose',0,'debug',0,'dualize',0);
res=optimize([LMI>=1*(1e-6),P>=0,lambda >=0],-trace(P),ops);
P = value(P);
lambda = value(lambda);
LMI = value(LMI);

feasible = res.problem == 0;
end

