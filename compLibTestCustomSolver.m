clear all
close all

% select a solver

% solver = 'customSolver';
% solver = 'Mosek';
solver = 'SeDuMi';
% solver = 'LMILab';

% select a list of problems

% problems = 'control_examples';
% problems = 'heat_flow_problems';
problems = 'second_order_problems';
load([problems,'.mat'])

% Warning:
% - Running any problem list with LMILab takes very long.
% - Running 'second_order_problems' or 'heat_flow_problems' also takes very
% long with the solvers 'Mosek', 'SeDuMi' and 'LMILab'.

% For the problems below took more 
% than 10^4 seconds.
% custom solver
% skip_examples = {};
% Mosek
% skip_examples = {'NN18','HF2D1','HF2D2','HF2D3','HF2D4','HF2D5','HF2D6','HF2D7','HF2D8','HF2D9','CM5','CM6','CBM'};
% SeDuMi
skip_examples = {'EB6','TL','NN18','HF2D1','HF2D2','HF2D3','HF2D4','HF2D5','HF2D6','HF2D7','HF2D8','HF2D9','CM4','CM5','CM6','CBM'};
% LMILab
% skip_examples = {'AC10','AC14','JE1','HF1','BDT2','CSE2','EB5','EB6','TL','CDP','NN18','HF2D1','HF2D2','HF2D3','HF2D4','HF2D5','HF2D6','HF2D7','HF2D8','HF2D9','CM2','CM3','CM4','CM5','CM6','DLR2','DLR3','ISS1','ISS2','CBM'};

for jjj = 1:length(model_list)
    name = model_list{jjj};
    
    addpath('./COMPlib_r1_1')

    %% Load a model from COMPlib
    [A,~,B,~,~,~,~,~,nx,~,nu] = COMPleib(name);
    % We consider only robust state feedback with input uncertainty. Thus,
    % we do not need most matrices in COMPleib.
    A = full(A);
    B1 = full(B);
    B2 = full(B);
    D1 = 0.25*eye(nu);
    C = zeros(nu,nx);
    D2 = zeros(nu,nu);
    nw = nu;
    nz = nu;

    %% Generate Control problem

    % For simplicity, we always study the robust LQR problem with
    Q = eye(nx);
    R = eye(nu);

    % We consider a robust state feedback synthesis problem with 25%
    % uncertainty in the input channels. To this end, we utilize diagonal
    % multipliers.

    s = nu;
    gamma = 1; % uncertainty gain

    M_tilde_gen = zeros(nw + nz,nw + nz,s);
    for iii = 1:s
        M_tilde_gen(iii,iii,iii) = -gamma^2;
        M_tilde_gen(nw+iii,nw+iii,iii) = 1;
    end

    %% Solve with custom solver
    disp(['Robust state feedback synthesis for problem ',name])
    disp(['Solver: ',solver])
    if any(strcmp(skip_examples,name))
        computation_times(jjj) = inf;
        feasibility_list(jjj) = false;
        achieved_cost(jjj) = inf;
        continue
    end
    tic
    if strcmp(solver,'customSolver')
        [lambda,P,LMI,feasible] = robustSynthesisCustomSolver(A,B1,B2,C,D1,D2,Q,R,M_tilde_gen);
    else
        [lambda,P,LMI,feasible] = robustSynthesisYalmip(A,B1,B2,C,D1,D2,Q,R,M_tilde_gen,solver);
    end
    computation_times(jjj) = toc;
    feasibility_list(jjj) = feasible;
    achieved_cost(jjj) = trace(P);
    disp(['Solution time: ',num2str(computation_times(jjj))])
    disp(['Feasible solution found: ',num2str(feasible)])
    disp(['Optimal value: ',num2str(achieved_cost(jjj))])
end

save([solver,'_',problems,'_results.mat'],'computation_times','feasibility_list','achieved_cost')