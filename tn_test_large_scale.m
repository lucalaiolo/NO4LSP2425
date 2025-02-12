%% Defining the variables for the test

clear
clc
close all

addpath('findiff_general')
addpath('solvers')
addpath('test_builds\exact_builds\')
addpath('test_builds\findiff_builds')

% Functions
[f1, gradf1, Hessf1] = exact_problem25();
[f2, gradf2, Hessf2] = exact_problem31();
[f3, gradf3, Hessf3] = exact_problem32();

% Forcing terms
fterms_lin = @(k, gradfk_norm) 0.5;
fterms_suplin = @(k, gradfk_norm) min(0.5, sqrt(gradfk_norm));
fterms_quad = @(k, gradfk_norm) min(0.5, gradfk_norm);

% Approximation choice
fterms_choice = input( ...
    'Choose a forcing term (1 for linear, 2 for superlinear, 3 for quadratic): ');

switch fterms_choice
    case 1
        fterms = fterms_lin;
    case 2
        fterms = fterms_suplin;
    case 3
        fterms = fterms_quad;
end

disp(' ');

disp('Choose type of finite difference approximation of the gradient');
disp('Possible choices:');
disp('fw -> forward difference approximation');
disp('c -> centered difference approximation');
disp('other -> no approximation');
FDgrad = input('Enter your choice: ', 's');

disp(' ');

disp('Choose type of finite difference approximation of the Hessian');
disp('Possible choices:');
disp('c -> centered difference approximation');
disp('Jfw -> forward difference approximation of Hessf as the Jacobian of the gradient');
disp('Jc -> centered difference approximation of Hessf as the Jacobian of the gradient');
disp('MF -> matrix-free implementation of the iterative method')
disp('other -> no approximation');
disp('!!! ATTENTION !!!')
disp('If FDgrad =  fw or FDgrad = c, the only possible choices are c and MF.');
disp('If a different option is selected, no approximation will be used.');
FDHess = input('Enter your choice: ', 's');
disp(' ');

seed = min([344611, 317663, 338344]);
rng(seed)
d = [3 4 5];
dim = 10.^d;
k = 2:2:8;
h_vec = 10.^-k;
% Value of the function at the computed solution
f_prob25 = zeros(11, length(dim), length(h_vec));
f_prob31 = zeros(11, length(dim), length(h_vec));
f_prob32 = zeros(11, length(dim), length(h_vec));

% Gradient norm computed at the computed solution
grad_norm_25 = zeros(11, length(dim), length(h_vec));
grad_norm_31 = zeros(11, length(dim), length(h_vec));
grad_norm_32 = zeros(11, length(dim), length(h_vec));

% Number of iterations for each problem
iterations_25 =  zeros(11, length(dim), length(h_vec));
iterations_31 = zeros(11, length(dim), length(h_vec));
iterations_32 = zeros(11, length(dim), length(h_vec));

% pcg iterations taken in each problem
pcg_it_25 = zeros(11, length(dim), length(h_vec));
pcg_it_31 = zeros(11, length(dim), length(h_vec));
pcg_it_32 = zeros(11, length(dim), length(h_vec));

% Execution times
time_25 = zeros(11, length(dim), length(h_vec));
time_31 = zeros(11, length(dim), length(h_vec));
time_32 = zeros(11, length(dim), length(h_vec));

% Exit flags that tell us if the method failed because it could not find a
% suitable steplength for the backtracking procedure
bcktrck_fail_25 = zeros(11, length(dim), length(h_vec));
bcktrck_fail_31 = zeros(11, length(dim), length(h_vec));
bcktrck_fail_32 = zeros(11, length(dim), length(h_vec));

% Distances from the suggested initial point
dist_25 = zeros(10, length(dim));
dist_31 = zeros(10, length(dim));
dist_32 = zeros(10, length(dim));

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
    
    dist_25(:, i) = vecnorm(x0_ros - x0_ros_rand)';
    dist_31(:, i) = vecnorm(x0_broyden - x0_broyden_rand);
    dist_32(:, i) = dist_31(:, i);

    % Solve the problems
    
    % Parameters
    kmax = 1e3;
    rho = 0.8;
    c1 = 1e-4;
    btmax = 50;
    pcg_maxit = 60;
    tolgrad = 1e-5;
    adapt = true; % To understand more about this, read findiff_grad (for example)
    
    for j = 1:length(h_vec)
        h = h_vec(j);
        
        tic;
        [~, f_prob25(1,i,j), grad_norm_25(1,i,j), bcktrck_fail_25(1,i,j), ...
                iterations_25(1,i,j),~,~,~, pcg_it_25(1,i,j)] = ...
            trunc_newton_general(x0_ros, f1, gradf1, Hessf1, kmax, ...
                tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, ...
                FDHess, h, adapt, false, ...
                @(gradf, x) findiff_Hess_problem25(gradf, x, h, adapt));
        time_25(1,i,j) = toc;
        %{
        tic;
        [~, f_prob31(1,i,j), grad_norm_31(1,i,j), bcktrck_fail_31(1,i,j), ...
                iterations_31(1,i,j),~,~,~, pcg_it_31(1,i,j)] = ...
            trunc_newton_general(x0_broyden, f2, gradf2, Hessf2, kmax, ...
                tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, ...
                FDHess, h, adapt, false, ...
                @(gradf, x) findiff_Hess_problem31_32(gradf, x, h, adapt));
        time_31(1,i,j) = toc;
        
        tic;
        [~, f_prob32(1,i,j), grad_norm_32(1,i,j), bcktrck_fail_32(1,i,j), ...
                iterations_32(1,i,j),~,~,~, pcg_it_32(1,i,j)] = ...
            trunc_newton_general(x0_broyden, f3, gradf3, Hessf3, kmax, ...
                tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, ...
                FDHess, h, adapt, true, ...
                @(gradf, x) findiff_Hess_problem31_32(gradf, x, h, adapt));
        time_32(1,i,j) = toc;
        %}
        for l = 1:10
            x0_1 = x0_ros_rand(:, l);
            x0_2 = x0_broyden_rand(:, l);
            x0_3 = x0_ext_broyden_rand(:, l);
            
            tic;
            [xk_1, f_prob25(l+1,i,j), grad_norm_25(l+1, i, j), bcktrck_fail_25(l+1, i, j), ...
                    iterations_25(l+1, i, j),~,~,~,pcg_it_25(l+1,i,j)] = ...
                trunc_newton_general(x0_1, f1, gradf1, Hessf1, kmax, ...
                    tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, ...
                    FDHess, h, adapt, false, ...
                    @(gradf, x) findiff_Hess_problem25(gradf, x, h, adapt));

            time_25(l+1, i, j) = toc;
            %{
            tic;
            [xk_2, f_prob31(l+1,i,j), grad_norm_31(l+1,i,j), bcktrck_fail_31(l+1,i,j), ... 
                    iterations_31(l+1,i,j),~,~,~,pcg_it_31(l+1,i,j)] = ...
                trunc_newton_general(x0_2, f2, gradf2, Hessf2, kmax, ...
                    tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, ...
                    FDHess, h, adapt, false, ...
                    @(gradf, x) findiff_Hess_problem31_32(gradf, x, h, adapt));

            time_31(l+1, i, j) = toc;
            disp(grad_norm_31(l+1,i,j));
            
            
            tic;
            [xk_3, f_prob32(l+1,i,j), grad_norm_32(l+1,i,j), bcktrck_fail_32(l+1,i,j), ...
                    iterations_32(l+1,i,j),~,~,~,pcg_it_32(l+1,i,j)] = ...
                trunc_newton_general(x0_3, f3, gradf3, Hessf3, kmax, ...
                    tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, ...
                    FDHess, h, adapt, true, ...
                    @(gradf, x) findiff_Hess_problem31_32(gradf, x, h, adapt));

            time_32(l+1, i, j) = toc;
            disp(grad_norm_32(l+1,i,j))
            %}
            
        end

    end

end


col1 = [10^3 * ones(44,1); 10^4 * ones(44,1); 10^5*ones(44,1)];
seq = zeros(44, 1);
seq(1:4) = 0;
seq(5:8) = 1;
seq(9:12) = 2;
seq(13:16)=3;
seq(17:20)=4;
seq(21:24)=5;
seq(25:28) = 6;
seq(29:32)=7;
seq(33:36) = 8;
seq(37:40) = 9;
seq(41:44) = 10;
col2 = [seq; seq; seq];

seq = repmat(h_vec', 11, 1);
col3 = [seq; seq; seq];

col4 = zeros(44*3, 1);
col5 = zeros(44*3, 1);
col6 = zeros(44*3, 1);
col7 = zeros(44*3, 1);
col8 = zeros(44*3, 1);
col9 = zeros(44*3, 1);
count = 1;
for dimensione = 1:3
    for test = 1:11
        for step = 1:4
            col4(count) = f_prob25(test, dimensione, step);
            col5(count) = grad_norm_25(test, dimensione, step);
            col6(count) = iterations_25(test, dimensione, step);
            col7(count) = pcg_it_25(test, dimensione, step);
            col8(count) = bcktrck_fail_25(test, dimensione, step);
            col9(count) = time_25(test, dimensione, step);
            count = count + 1;
        end
    end
end
T = table(col1, col2, col3, col4, col5, col6, col7, col8, col9);

table2latex(T(1:44,:), 'pag1');
table2latex(T(45:88,:), 'pag2');
table2latex(T(89:132,:), 'pag3');

%{
col1 = [10^3 * ones(11,1); 10^4 * ones(11,1); 10^5*ones(11,1)];
seq = 0:10;
col2 = [seq';seq';seq'];
col3 = 1*ones(33,1); % to change
col4 = [0;dist_31(:, 1);0; dist_31(:, 2);0; dist_31(:, 3)];
col5 = [f_prob31(:, 1); f_prob31(:, 2); f_prob31(:, 3)];
col6 = [grad_norm_31(:, 1); grad_norm_31(:, 2); grad_norm_31(:, 3)];
col7 = [iterations_31(:, 1); iterations_31(:, 2); iterations_31(:, 3)];
col8 = [pcg_it_31(:, 1); pcg_it_31(:, 2); pcg_it_31(:, 3)];
col9 = [bcktrck_fail_31(:, 1); bcktrck_fail_31(:, 2); bcktrck_fail_31(:, 3)];
col10 = [time_31(:, 1); time_31(:, 2); time_31(:, 3)];

T = table(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10);

table2latex(T);


col1 = [10^3 * ones(11,1); 10^4 * ones(11,1); 10^5*ones(11,1)];
seq = 0:10;
col2 = [seq';seq';seq'];
col3 = 2*ones(33,1); % to change
col4 = [0;dist_32(:, 1);0; dist_32(:, 2);0; dist_32(:, 3)];
col5 = [f_prob32(:, 1); f_prob32(:, 2); f_prob32(:, 3)];
col6 = [grad_norm_32(:, 1); grad_norm_32(:, 2); grad_norm_32(:, 3)];
col7 = [iterations_32(:, 1); iterations_32(:, 2); iterations_32(:, 3)];
col8 = [pcg_it_32(:, 1); pcg_it_32(:, 2); pcg_it_32(:, 3)];
col9 = [bcktrck_fail_32(:, 1); bcktrck_fail_32(:, 2); bcktrck_fail_32(:, 3)];
col10 = [time_32(:, 1); time_32(:, 2); time_32(:, 3)];

T = table(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10);

table2latex(T);

%}