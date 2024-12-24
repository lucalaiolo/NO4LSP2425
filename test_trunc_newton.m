%% Defining the variables for the test

clear
clc
close all

% Functions
[f1, gradf1, Hessf1] = exact_rosenbrock();
[f2, gradf2, Hessf2] = exact_ext_rosenbrock();
[f3, gradf3, Hessf3] = exact_banded_trigon();

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
disp('In the case of a large scale problem, it is suggested to use MF.')
FDHess = input('Enter your choice: ', 's');
disp(' ');

seed = min([9876, 1234, 2345]); % to do

d = 3:5;
dim = 10.^d;
k = 6:2:12;
h_vec = 10.^-k;

% Number of iterations for each example
iterations_ros =  zeros(10, length(dim), length(h_vec));
iterations_ext_ros = zeros(10, length(dim), length(h_vec));
iterations_trig = zeros(10, length(dim), length(h_vec));

% Number of successes for each problem and for each dimension
success_ros = zeros(1, length(dim), length(h_vec));
success_ext_ros = zeros(1, length(dim), length(h_vec));
success_trig = zeros(1, length(dim), length(h_vec));

% Execution times
time_ros = zeros(10, length(dim), length(h_vec));
time_ext_ros = zeros(10, length(dim), length(h_vec));
time_trig = zeros(10, length(dim), length(h_vec));

% Orders of convergence 

% to do

% Distance from the solution
dist_1 = zeros(10, length(dim), length(h_vec));
dist_2 = zeros(10, length(dim), length(h_vec));
% dist_3 = to do


%% Running the tests

for i = 1:length(dim)
    n = dim(i);

    % Rosenbrock suggested starting point
    x0_rosenbrock = ones(n, 1);
    x0_rosenbrock(1:2:end) = -1.2;

    % Extended Rosenbrock suggested starting point
    x0_ext_rosenbrock = x0_rosenbrock; % same value

    % Banded trigonometric problem suggested starting point
    x0_trig = ones(n, 1);

    % Perturbed initial conditions
    x0_rosenbrock_rand = x0_rosenbrock + 2 .* rand(n, 10) - 1;
    x0_ext_rosenbrock_rand = x0_ext_rosenbrock + 2 .* rand(n, 10) - 1;
    x0_trig_rand = x0_trig + 2 .* rand(n, 10) - 1;

    % Solve the problems
    
    % Parameters
    kmax = 10000;
    rho = 0.8;
    c1 = 1e-4;
    btmax = 75;
    pcg_maxit = 50;
    tolgrad = 1e-7;
    adapt = true; % To understand more about this, read findiff_grad (for example)
    
    for j = 1:length(h_vec)
        h = h_vec(j);
        for l = 1:10
            x0_1 = x0_rosenbrock_rand(:, l);
            x0_2 = x0_ext_rosenbrock_rand(:, l);
            x0_3 = x0_trig_rand(:, l);
            %{
            tic;
            [xk_1, fk_1, gradfk_norm_1, k_1, xseq_1, btseq_1, flag_1, pcg_iterseq_1] = ...
                trunc_newton_general(x0_1, f1, gradf1, Hessf1, kmax, ...
                tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, ...
                FDHess, h, adapt, true, @(gradf, x) findiff_Hess_tridiagonal(gradf, x, h));
            time_ros(l, i, j) = toc;
            %}
            
            tic;
            [xk_2, fk_2, gradfk_norm_2, k_2, xseq_2, btseq_2, flag_2, pcg_iterseq_2] = ...
                trunc_newton_general(x0_2, f2, gradf2, Hessf2, kmax, ...
                tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, ...
                FDHess, h, adapt, true, @(gradf, x) findiff_Hess_ext_rosenbrock(gradf, x, h));
            time_ext_ros(l, i, j) = toc;
            dist_2(l, i, j) = norm(xk_2 - ones(n, 1));
            
            %{ 
            tic;
            [xk_3, fk_3, gradfk_norm_3, k_3, xseq_3, btseq_3, flag_3, pcg_iterseq_3] = ...
                trunc_newton_general(x0_3, f3, gradf3, Hessf3, kmax, ...
                tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, ...
                FDHess, h, adapt, true, @(gradf, x) findiff_Hess_diagonal(gradf, x, h));
            time_trig(l, i, j) = toc;
            %}
        end

    end

end