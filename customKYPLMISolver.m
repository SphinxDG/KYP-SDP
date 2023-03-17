function [lambda,lambda0,Kp,Pp,Q,S,R] = customKYPLMISolver(c,Sigma,A,B,Q_gen,S_gen,R_gen,N_gen,lambda)

[n,m] = size(B);
s = size(Q_gen,3); % the last multiplier is the feasibility multiplier
n_lam_active = s; % in the beginning all multipliers are optimized. Later, the feasibility multiplier is no longer used.
Q = sum(pagemtimes(Q_gen,permute([1;lambda],[2,3,1])),3);
S = sum(pagemtimes(S_gen,permute([1;lambda],[2,3,1])),3);
R = sum(pagemtimes(R_gen,permute([1;lambda],[2,3,1])),3);
Q_gen(:,:,s+1) = -eye(n);
R_gen(:,:,s+1) = -eye(m);
S_gen(:,:,s+1) = zeros(n,m);
N_gen(:,:,s+1) = zeros(size(N_gen,1));
lambda(s,1) = max(max(real(eig([Q,S;S',R])))*1.1,max(real(eig([Q,S;S',R])))+0.1);
if lambda(s) > 0
    infeasible = true;
else
    infeasible = false;
end

keep_optimizing = true;

t = inf;
[feasible,v,Pp,Pm,Kp,Km,Q,S,R,N] = evaluateRiccati(c,Sigma,lambda,Q_gen,S_gen,R_gen,N_gen,A,B,t,infeasible);

if ~feasible
    keep_optimizing = false;
end

dlogdet = @(mat,dmat,ind) - trace(mat\dmat(:,:,ind));
d2logdet = @ (mat,dmat,d2mat,ind1,ind2) trace((mat\dmat(:,:,ind1))*(mat\dmat(:,:,ind2))) - trace(mat\d2mat(:,:,ind1,ind2));
num_Riccati = 1;
while keep_optimizing
    
    dv = zeros(n_lam_active,1);
    d2v = zeros(n_lam_active,n_lam_active);
    [dPp,d2Pp] = differentiateThroughRiccati(Kp,A,B,R,S_gen,Q_gen,R_gen,n_lam_active);
    [dPm,d2Pm] = differentiateThroughRiccati(Km,A,B,R,S_gen,Q_gen,R_gen,n_lam_active);
    
    % gradient computation
    for i2 = 1:n_lam_active
        dv(i2) = dlogdet(Pp,dPp,i2) + dlogdet(Pp-Pm,dPp-dPm,i2) + dlogdet(N,N_gen(:,:,2:end),i2) + dlogdet(-R , -R_gen(:,:,2:end),i2);
    end
    
    % hessian computation
    for i2 = 1:n_lam_active
        for i3 = 1:i2
            d2v(i2,i3) = d2logdet(Pp,dPp,d2Pp,i2,i3) + d2logdet(Pp-Pm,dPp-dPm,d2Pp-d2Pm,i2,i3) + trace((R\R_gen(:,:,1+i2))*(R\R_gen(:,:,1+i3))) + trace((N\N_gen(:,:,1+i2))*(N\N_gen(:,:,1+i3)));
            d2v(i3,i2) = d2v(i2,i3);
        end
    end
    if infeasible
        if t == inf
            t = -dv(end);
        end
        dv(end) = dv(end) + t;
    else
        dv_add = -permute(sum(Sigma.*dPp,[1,2]),[3,1,2]);
        dv_add(1:n_lam_active) = dv_add(1:n_lam_active) + c;
        d2v_add = -permute(sum(Sigma.*d2Pp,[1,2]),[3,4,1,2]);
        if t == inf
            t = norm(dv)/norm(dv_add);
        end
        dv = dv + t*dv_add(1:n_lam_active);
        d2v = d2v + t*d2v_add(1:n_lam_active,1:n_lam_active);
    end
    
    % step size control
    Newton_dec = sqrt(dv(1:n_lam_active)'*(d2v(1:n_lam_active,1:n_lam_active)\dv(1:n_lam_active)));
    if Newton_dec > 0.25
        accuracy = 0;
        alpha = 1;
    else
        accuracy = Newton_dec^2 + (n+n+m+size(N_gen,1))/t;
        alpha = 1;
    end
    line_search_success = false;
    while true
        step = -alpha*(d2v\dv);
        lambda_alpha(1:n_lam_active,1) = lambda(1:n_lam_active) + step;
        [feasible,v_alpha,Pp_alpha,Pm_alpha,Kp_alpha,Km_alpha,Q_alpha,S_alpha,R_alpha,N_alpha] = evaluateRiccati(c,Sigma,lambda_alpha,Q_gen,S_gen,R_gen,N_gen,A,B,t,infeasible);
        num_Riccati = num_Riccati + ~isempty(Pp);
        if feasible && (v_alpha <= v + (step'*dv(1:n_lam_active))*(alpha/4))
            v = v_alpha; Pp = Pp_alpha; Pm = Pm_alpha; Kp = Kp_alpha; Km = Km_alpha; Q = Q_alpha; S = S_alpha; R = R_alpha; N = N_alpha;
            lambda = lambda_alpha;
            line_search_success = true;
            break
        elseif abs(step'*dv(1:n_lam_active)) <= 100*eps
            break
        else
            alpha = min(alpha/2,sqrt(alpha)/sqrt(1+Newton_dec));
        end
    end
    if infeasible && lambda(end) <= 0
        infeasible = false;
        n_lam_active = s-1;
        t = inf;
        [~,v,Pp,Pm,Kp,Km,Q,S,R,N] = evaluateRiccati(c,Sigma,lambda,Q_gen,S_gen,R_gen,N_gen,A,B,t,infeasible);
        num_Riccati = num_Riccati + 1;
    elseif (lambda(end) - accuracy >= 0) && accuracy > 0
            keep_optimizing = false;
    elseif infeasible && (Newton_dec^2 <= 0.25*accuracy || ~line_search_success || dv(end) < 0)
        if t <= 1e10
            v = v - t*lambda(end);
            t = max(10*t,t-dv(end));
            v = v + t*lambda(end);
        elseif ~line_search_success
            keep_optimizing = false;
        end
    elseif ~infeasible && (Newton_dec^2 <= 0.25*accuracy || ~line_search_success)
        if accuracy <= (1e-6)*abs(v)
            keep_optimizing = false;
        elseif t <= 1e10
            v = v - t*(c'*lambda(1:end-1) - trace(Sigma*Pp));
            t = 10*t;
            v = v + t*(c'*lambda(1:end-1) - trace(Sigma*Pp));
        elseif ~line_search_success
            keep_optimizing = false;
        end
    end
end

%% Retrieve multipliers

lambda0 = lambda(end);
lambda = lambda(1:end-1);

num_Riccati
end

function [feasible,v,Pp,Pm,Kp,Km,Q,S,R,N] = evaluateRiccati(c,Sigma,lambda,Q_gen,S_gen,R_gen,N_gen,A,B,t,infeasible)
feasible = true;
v = [];
Pp = [];
Kp = [];
Pm = [];
Km = [];
Q = [];
S = [];
R = [];
n = size(A,1);
s = size(lambda,1);
N = sum(pagemtimes(N_gen,permute([1;lambda],[2,3,1])),3);
[cholN,FLAG] = chol(N);
if FLAG ~= 0
    feasible = false;
    return
end
if infeasible
    R = sum(pagemtimes(R_gen,permute([1;lambda],[2,3,1])),3);
else
    R = sum(pagemtimes(R_gen(:,:,1:s),permute([1;lambda(1:s-1)],[2,3,1])),3);
end
[cholR,FLAG] = chol(-R);
if FLAG ~= 0
    feasible = false;
    return
end
if infeasible
    Q = sum(pagemtimes(Q_gen,permute([1;lambda],[2,3,1])),3);
    S = sum(pagemtimes(S_gen,permute([1;lambda],[2,3,1])),3);
else
    Q = sum(pagemtimes(Q_gen(:,:,1:s),permute([1;lambda(1:s-1)],[2,3,1])),3);
    S = sum(pagemtimes(S_gen(:,:,1:s),permute([1;lambda(1:s-1)],[2,3,1])),3);
end

Pp = careSchur(A,B,Q,S,R);
[cholP,FLAG] = chol(Pp);
if FLAG ~= 0 || isempty(Pp)
    feasible = false;
    return
end
Kp = R\(S+Pp*B)';
DPinv = lyap((A-B*Kp),(B/R)*B');
[cholDPinv,FLAG] = chol(DPinv);
if FLAG ~= 0
    feasible = false;
    return
end
cholDP = inv(cholDPinv);
Pm = Pp-cholDP*cholDP';
Km = R\(S+Pm*B)';
b = @(cholmat) -2*sum(log(diag(cholmat)));
v = b(cholP) + b(cholDP) + b(cholN) + b(cholR);
if t == inf
    v = inf;
elseif infeasible
    v = v + t*lambda(end);
else
    v = v + t*(c'*lambda(1:end-1) - trace(Sigma*Pp));
end
end

function [dP,d2P] = differentiateThroughRiccati(K,A,B,R,S_gen,Q_gen,R_gen,n_lam_active)

[n,m,s] = size(S_gen);
s = s-1;

dP = zeros(n,n,s);
dK = zeros(m,n,s);
d2P = zeros(n,n,s,s);

Aanti = (A-B*K)';

for i2 = 1:n_lam_active
    dP(:,:,i2) = lyap(Aanti,[eye(n);-K]'*[Q_gen(:,:,1+i2),S_gen(:,:,1+i2);S_gen(:,:,1+i2)',R_gen(:,:,1+i2)]*[eye(n);-K]);
    dK(:,:,i2) = R\(B'*dP(:,:,i2) + S_gen(:,:,1+i2)' - R_gen(:,:,1+i2)*K);
end

% hessian computation
for i2 = 1:n_lam_active
    for i3 = 1:i2
        dKRdK = dK(:,:,i2)'*R*dK(:,:,i3);
        d2P(:,:,i2,i3) = lyap(Aanti,-dKRdK -dKRdK');
        d2P(:,:,i3,i2) = d2P(:,:,i2,i3);
    end
end
    
end

function P = careSchur(A,B,Q,S,R)

H = [A - (B/R)*S',-(B/R)*B';
     (S/R)*S'-Q, S*(R\B') - A'];

[T1,M1] = schur(H);

if sum(diag(M1)>0) == sum(diag(M1)<0)
    [T,~] = ordschur(T1,M1,'rhp');
    
    n = size(A,1);
    V11 = T(1:n,1:n);
    V21 = T(n+(1:n),1:n);
    
    P = V21/V11;
    P = 0.5*(P + P');
else
    P = [];
end
end