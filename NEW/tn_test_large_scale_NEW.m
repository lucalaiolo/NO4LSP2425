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
k = 2:2:12;
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
    adapt = false; % To understand more about this, read findiff_grad (for example)
    
    for j = 1:length(h_vec)
        h = h_vec(j);
        
        tic;
        [~, f_prob25(1,i,j), grad_norm_25(1,i,j), bcktrck_fail_25(1,i,j), ...
                iterations_25(1,i,j),~,~,~, pcg_it_25(1,i,j)] = ...
            trunc_newton_general_NEW(x0_ros, f1, gradf1, Hessf1, kmax, ...
                tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, ...
                FDHess, h, adapt, true, @(x) findiff_grad_problem25(f1, x, h, adapt), ...
                true, @(gradf, x) findiff_Hess_problem25(gradf, x, h, adapt));
        time_25(1,i,j) = toc;
        %{
        tic;
        [~, f_prob31(1,i,j), grad_norm_31(1,i,j), bcktrck_fail_31(1,i,j), ...
                iterations_31(1,i,j),~,~,~, pcg_it_31(1,i,j)] = ...
            trunc_newton_general_NEW(x0_broyden, f2, gradf2, Hessf2, kmax, ...
                tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, ...
                FDHess, h, adapt, true, @(x) findiff_grad_problem31(f2, x, h, adapt), ...
                true, @(gradf, x) findiff_Hess_problem31_32(gradf, x, h, adapt));
        time_31(1,i,j) = toc;
        
        tic;
        [~, f_prob32(1,i,j), grad_norm_32(1,i,j), bcktrck_fail_32(1,i,j), ...
                iterations_32(1,i,j),~,~,~, pcg_it_32(1,i,j)] = ...
            trunc_newton_general_NEW(x0_broyden, f3, gradf3, Hessf3, kmax, ...
                tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, ...
                FDHess, h, adapt, true, @(x) findiff_grad_problem32(f3, x, h, adapt), ...
                false, @(gradf, x) findiff_Hess_problem31_32(gradf, x, h, adapt));
        time_32(1,i,j) = toc;
        %}
        for l = 1:10
            x0_1 = x0_ros_rand(:, l);
            x0_2 = x0_broyden_rand(:, l);
            x0_3 = x0_ext_broyden_rand(:, l);
            
            tic;
            [xk_1, f_prob25(l+1,i,j), grad_norm_25(l+1, i, j), bcktrck_fail_25(l+1, i, j), ...
                    iterations_25(l+1, i, j),~,~,~,pcg_it_25(l+1,i,j)] = ...
                trunc_newton_general_NEW(x0_1, f1, gradf1, Hessf1, kmax, ...
                    tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, ...
                    FDHess, h, adapt, true, @(x) findiff_grad_problem25(f1, x, h, adapt), ...
                    true, @(gradf, x) findiff_Hess_problem25(gradf, x, h, adapt));
            
            time_25(l+1, i, j) = toc;
            disp(grad_norm_25(l+1, i, j));
            %{
            tic;
            [xk_2, f_prob31(l+1,i,j), grad_norm_31(l+1,i,j), bcktrck_fail_31(l+1,i,j), ... 
                    iterations_31(l+1,i,j),~,~,~,pcg_it_31(l+1,i,j)] = ...
                trunc_newton_general_NEW(x0_2, f2, gradf2, Hessf2, kmax, ...
                    tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, ...
                    FDHess, h, adapt, true, @(x) findiff_grad_problem31(f2, x, h, adapt), ...
                    true, @(gradf, x) findiff_Hess_problem31_32(gradf, x, h, adapt));

            time_31(l+1, i, j) = toc;
            disp(grad_norm_31(l+1,i,j));
            
            
            tic;
            [xk_3, f_prob32(l+1,i,j), grad_norm_32(l+1,i,j), bcktrck_fail_32(l+1,i,j), ...
                    iterations_32(l+1,i,j),~,~,~,pcg_it_32(l+1,i,j)] = ...
                trunc_newton_general_NEW(x0_3, f3, gradf3, Hessf3, kmax, ...
                    tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, ...
                    FDHess, h, adapt, true, @(x) findiff_grad_problem32(f3, x, h, adapt), ...
                    false, @(gradf, x) findiff_Hess_problem31_32(gradf, x, h, adapt));

            time_32(l+1, i, j) = toc;
            disp(grad_norm_32(l+1,i,j))
            %}
            
        end

    end

end
