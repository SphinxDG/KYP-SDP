clear all
close all

%% Computation time table

solvers = {'customSolver','Mosek','SeDuMi','LMILab'};
problem_lists = {'control_examples','heat_flow_problems','second_order_problems'};

time_struct = struct();
time_struct.Problem = [];

for ii = 1:lenth(solvers)
    solver = solvers{ii};
    time_struct.(solver) = [];
    for jj = 1:lenth(problem_lists)
        problems = problem_lists{jj};
        load([solver,'_',problems,'_results.mat']);
        load([problems,'.mat'])
        time_struct.Problem = {time_struct.Problem,model_list};
        time_struct.(solver) = [time_struct.(solver);computation_times];
    end
end