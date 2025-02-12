%% Defining the variables for the test

clear
clc
close all

addpath('solvers')
addpath('test_builds\exact_builds\')

% Functions
[f1, gradf1] = exact_problem25();
[f2, gradf2] = exact_problem31();
[f3, gradf3] = exact_problem32();

seed = min([344611, 317663, 338344]);
rng(seed)
dim = [10 25 50]';

% Value of the function at the computed solution
f_prob25 = zeros(10, length(dim));
f_prob31 = zeros(10, length(dim));
f_prob32 = zeros(10, length(dim));

% Gradient norm computed at the computed solution
grad_norm_25 = zeros(10, length(dim));
grad_norm_31 = zeros(10, length(dim));
grad_norm_32 = zeros(10, length(dim));

% Number of iterations for each problem
iterations_25 =  zeros(10, length(dim));
iterations_31 = zeros(10, length(dim));
iterations_32 = zeros(10, length(dim));

% Number of function evaluations for each problem
n_fun_evals_25 = zeros(10, length(dim));
n_fun_evals_31 = zeros(10, length(dim));
n_fun_evals_32 = zeros(10, length(dim));

% Execution times
time_25 = zeros(10, length(dim));
time_31 = zeros(10, length(dim));
time_32 = zeros(10, length(dim));

% Number of transformations for each problem
transf_25 = zeros(10, length(dim), 4);
transf_31 = zeros(10, length(dim), 4);
transf_32 = zeros(10, length(dim), 4);

%% Running the tests

for i = 1:length(dim)
    n = dim(i);

    % Problem 25 suggested starting point
    x0_ros = ones(n, 1);
    x0_ros(1:2:end) = -1.2;

    % Problem 31 suggested starting point
    x0_broyden = -ones(n, 1);

    % Problem 32 suggested starting point
    x0_ext_broyden = x0_broyden;

    % Perturbed initial conditions
    x0_ros_rand = x0_ros + 2 .* rand(n, 10) - 1;
    x0_broyden_rand = x0_broyden + 2 .* rand(n, 10) - 1;
    x0_ext_broyden_rand = x0_broyden_rand;

    % Solve the problems
    
    % Parameters
    rho = 1;
    chi = 2;
    gamma = 0.5;
    sigma = 0.5;
    tolFun = 1e-12;
    tolX = 1e-12;
    MaxIter = 12000 * n;
    MaxFunEvals = 16000 * n;
    c = 0.05;
    for l = 1:10
        x0_1 = x0_ros_rand(:, l);
        x0_2 = x0_broyden_rand(:, l);
        x0_3 = x0_ext_broyden_rand(:, l);
        disp(n);
        %{
        tic;
        
        [sol_1, iterations_25(l, i), n_fun_evals_25(l, i), f_seq_1, ...
                transformation_1, best_seq_1] = ...
            nelder_mead(f1, x0_1, tolFun, tolX, MaxIter, MaxFunEvals, ...
                c, rho, chi, gamma, sigma);

        time_25(l, i) = toc;
        f_prob25(l, i) = f_seq_1(end);
        grad_norm_25(l, i) = norm(gradf1(sol_1));
        disp(f_seq_1(end));
        disp(grad_norm_25(l,i));
        for k=1:length(transformation_1)
            switch transformation_1{k}
                case 'expansion'
                    transf_25(l, i, 1) = transf_25(l, i, 1) + 1;
                case 'reflection'
                    transf_25(l, i, 2) = transf_25(l, i, 2) + 1;
                case 'contraction'
                    transf_25(l, i, 3) = transf_25(l, i, 3) + 1;
                case 'shrinkage'
                    transf_25(l, i, 4) = transf_25(l, i, 4) + 1;
            end
        end
        %}
        %{
        tic;
        
        [sol_2, iterations_31(l, i), n_fun_evals_31(l, i), f_seq_2, ...
                transformation_2, best_seq_2] = ...
            nelder_mead(f2, x0_2, tolFun, tolX, MaxIter, MaxFunEvals, ...
                c, rho, chi, gamma, sigma);

        time_31(l, i) = toc;
        f_prob31(l, i) = f_seq_2(end);
        grad_norm_31(l, i) = norm(gradf2(sol_2));
        disp(f_seq_2(end));
        disp(grad_norm_31(l,i));
        for k=1:length(transformation_2)
            switch transformation_2{k}
                case 'expansion'
                    transf_31(l, i, 1) = transf_31(l, i, 1) + 1;
                case 'reflection'
                    transf_31(l, i, 2) = transf_31(l, i, 2) + 1;
                case 'contraction'
                    transf_31(l, i, 3) = transf_31(l, i, 3) + 1;
                case 'shrinkage'
                    transf_31(l, i, 4) = transf_31(l, i, 4) + 1;
            end
        end
        %}
        
        tic;
        
        [sol_3, iterations_32(l, i), n_fun_evals_32(l, i), f_seq_3, ...
                transformation_3, best_seq_3] = ...
            nelder_mead(f3, x0_3, tolFun, tolX, MaxIter, MaxFunEvals, ...
                c, rho, chi, gamma, sigma);

        time_32(l, i) = toc;
        f_prob32(l, i) = f_seq_3(end);
        grad_norm_32(l, i) = norm(gradf3(sol_3));
        disp(f_seq_3(end));
        disp(grad_norm_32(l,i));
        for k=1:length(transformation_3)
            switch transformation_3{k}
                case 'expansion'
                    transf_32(l, i, 1) = transf_32(l, i, 1) + 1;
                case 'reflection'
                    transf_32(l, i, 2) = transf_32(l, i, 2) + 1;
                case 'contraction'
                    transf_32(l, i, 3) = transf_32(l, i, 3) + 1;
                case 'shrinkage'
                    transf_32(l, i, 4) = transf_32(l, i, 4) + 1;
            end
        end
        
    end


end
%{
col1 = [10*ones(10,1); 25 * ones(10,1); 50 * ones(10,1)];
seq = 1:10;
col2 = [seq'; seq'; seq'];
col3 = [f_prob25(1:10, 1); f_prob25(1:10, 2); f_prob25(1:10, 3)];
col4 = [grad_norm_25(1:10, 1); grad_norm_25(1:10, 2); grad_norm_25(1:10, 3)];
col5 = [iterations_25(1:10, 1); iterations_25(1:10, 2); iterations_25(1:10, 3)];
col6 = [n_fun_evals_25(1:10, 1); n_fun_evals_25(1:10, 2); n_fun_evals_25(1:10, 3)];
col7 = [transf_25(1:10, 1, 1); transf_25(1:10, 2, 1); transf_25(1:10, 3, 1)];
col8 = [transf_25(1:10, 1, 2); transf_25(1:10, 2, 2); transf_25(1:10, 3, 2)];
col9 = [transf_25(1:10, 1, 3); transf_25(1:10, 2, 3); transf_25(1:10, 3, 3)];
col10 = [transf_25(1:10, 1, 4); transf_25(1:10, 2, 4); transf_25(1:10, 3, 4)];
col11 = [time_25(1:10, 1); time_25(1:10, 2); time_25(1:10, 3)];
T = table(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11);
table2latex(T);
%}
%{
col1 = [10*ones(10,1); 25 * ones(10,1); 50 * ones(10,1)];
seq = 1:10;
col2 = [seq'; seq'; seq'];
col3 = [f_prob31(1:10, 1); f_prob31(1:10, 2); f_prob31(1:10, 3)];
col4 = [grad_norm_31(1:10, 1); grad_norm_31(1:10, 2); grad_norm_31(1:10, 3)];
col5 = [iterations_31(1:10, 1); iterations_31(1:10, 2); iterations_31(1:10, 3)];
col6 = [n_fun_evals_31(1:10, 1); n_fun_evals_31(1:10, 2); n_fun_evals_31(1:10, 3)];
col7 = [transf_31(1:10, 1, 1); transf_31(1:10, 2, 1); transf_31(1:10, 3, 1)];
col8 = [transf_31(1:10, 1, 2); transf_31(1:10, 2, 2); transf_31(1:10, 3, 2)];
col9 = [transf_31(1:10, 1, 3); transf_31(1:10, 2, 3); transf_31(1:10, 3, 3)];
col10 = [transf_31(1:10, 1, 4); transf_31(1:10, 2, 4); transf_31(1:10, 3, 4)];
col11 = [time_31(1:10, 1); time_31(1:10, 2); time_31(1:10, 3)];
T = table(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11);
table2latex(T);
%}

col1 = [10*ones(10,1); 25 * ones(10,1); 50 * ones(10,1)];
seq = 1:10;
col2 = [seq'; seq'; seq'];
col3 = [f_prob32(1:10, 1); f_prob32(1:10, 2); f_prob32(1:10, 3)];
col4 = [grad_norm_32(1:10, 1); grad_norm_32(1:10, 2); grad_norm_32(1:10, 3)];
col5 = [iterations_32(1:10, 1); iterations_32(1:10, 2); iterations_32(1:10, 3)];
col6 = [n_fun_evals_32(1:10, 1); n_fun_evals_32(1:10, 2); n_fun_evals_32(1:10, 3)];
col7 = [transf_32(1:10, 1, 1); transf_32(1:10, 2, 1); transf_32(1:10, 3, 1)];
col8 = [transf_32(1:10, 1, 2); transf_32(1:10, 2, 2); transf_32(1:10, 3, 2)];
col9 = [transf_32(1:10, 1, 3); transf_32(1:10, 2, 3); transf_32(1:10, 3, 3)];
col10 = [transf_32(1:10, 1, 4); transf_32(1:10, 2, 4); transf_32(1:10, 3, 4)];
col11 = [time_32(1:10, 1); time_32(1:10, 2); time_32(1:10, 3)];
T = table(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11);
table2latex(T);
