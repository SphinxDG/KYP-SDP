clear all
close all

%% Computation time table

solvers = {'customSolver','Mosek','SeDuMi','LMILab'};
problem_lists = {'control_examples','heat_flow_problems','second_order_problems'};

time_struct = struct();
time_struct.Problem = [];

for ii = 1:length(solvers)
    solver = solvers{ii};
    time_struct.(solver) = [];
    for jj = 1:length(problem_lists)
        problems = problem_lists{jj};
        load([solver,'_',problems,'_results.mat']);
        time_struct.(solver) = [time_struct.(solver);computation_times'];
    end
end

for jj = 1:length(problem_lists)
    problems = problem_lists{jj};
    load([problems,'.mat'])
    time_struct.Problem = [time_struct.Problem;model_list'];
end

computationTimes_table = struct2table(time_struct)
% table2latex(computationTimes_table,'computationTimes.tex')


%% Optimal value table

solvers = {'customSolver','Mosek','SeDuMi','LMILab'};
problem_lists = {'control_examples','heat_flow_problems','second_order_problems'};

value_struct = struct();
value_struct.Problem = [];

for ii = 1:length(solvers)
    solver = solvers{ii};
    value_struct.(solver) = [];
    for jj = 1:length(problem_lists)
        problems = problem_lists{jj};
        load([solver,'_',problems,'_results.mat']);
        value_struct.(solver) = [value_struct.(solver);achieved_cost'];
    end
end

for jj = 1:length(problem_lists)
    problems = problem_lists{jj};
    load([problems,'.mat'])
    value_struct.Problem = [value_struct.Problem;model_list'];
end

value_table = struct2table(value_struct)
% table2latex(value_table,'optimalValues.tex')