%% Defining the variables for the test

clear
clc
close all

% f, gradf, Hessf

f = @(x) 100 * (x(2,:) - x(1,:).^2).^2 + (1 - x(1,:)).^2;

gradf = @(x) [400 * x(1,:).^3 - 400 * x(1,:) .* x(2,:) + 2 * x(1,:) - 2; ...
    200* (x(2,:) - x(1,:).^2)];

Hessf = @(x) [1200 * x(1)^2 - 400 * x(2) + 2, -400 * x(1); -400 * x(1), 200];

% Forcing terms

fterms_lin = @(k, gradfk_norm) 0.5;
fterms_suplin = @(k, gradfk_norm) min(0.5, sqrt(gradfk_norm));
fterms_quad = @(k, gradfk_norm) min(0.5, gradfk_norm);

% Starting points for the test

x0_1 = [1.2; 1.2];
x0_2 = [-1.2; 1];

% Default choices for the truncated Newton method
kmax = 10000;
tolgrad = 1e-8;
pcg_maxit = 500;
c1 = 1e-4;
rho = 0.8;
btmax = 75;
h = sqrt(eps);
adapt = true; % To understand more about this, read findiff_grad (for example)

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

%% Run the tests
disp('Test 1: x0 = [1.2, 1.2]')

disp('**** TRUNCATED NEWTON (BACKTRACK + FLAG): START *****')
[xk_tn_1, fk_tn_1, gradfk_norm_tn_1, k_tn_1, xseq_tn_1, btseq_tn_1, flag_tn_1, ...
    pcg_iterseq_tn_1] = trunc_newton_general(x0_1, f, gradf, Hessf, kmax, ...
        tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, FDHess, h, adapt);
disp('**** TRUNCATED NEWTON: FINISHED *****')
disp('**** TRUNCATED NEWTON: RESULTS *****')
disp('************************************')
disp(['xk: ', mat2str(xk_tn_1), ' (actual solution: [1; 1])'])
disp(['f(xk): ', num2str(fk_tn_1), ' (actual minimum: 0)'])
disp('FLAG CONTENT:')
disp(flag_tn_1)
disp('************************************')

disp(' ')

disp('Test 2: x0 = [-1.2; 1]')
disp('**** TRUNCATED NEWTON (BACKTRACK + FLAG): START *****')
[xk_tn_2, fk_tn_2, gradfk_norm_tn_2, k_tn_2, xseq_tn_2, btseq_tn_2, flag_tn_2, ...
    pcg_iterseq_tn_2] = trunc_newton_general(x0_2, f, gradf, Hessf, kmax, ...
        tolgrad, pcg_maxit, fterms, c1, rho, btmax, FDgrad, FDHess, h, adapt);
disp('**** TRUNCATED NEWTON: FINISHED *****')
disp('**** TRUNCATED NEWTON: RESULTS *****')
disp('************************************')
disp(['xk: ', mat2str(xk_tn_2), ' (actual solution: [1; 1])'])
disp(['f(xk): ', num2str(fk_tn_2), ' (actual minimum: 0)'])
disp('FLAG CONTENT:')
disp(flag_tn_2)
disp('************************************')

%% Plots

f_meshgrid = @(X,Y)reshape(f([X(:),Y(:)]'),size(X));
% Creation of the meshgrid for the contour-plot
[X, Y] = meshgrid(linspace(-6, 6, 500), linspace(-6, 6, 500));
% Computation of the values of f for each point of the mesh
Z = f_meshgrid(X, Y);
fig1 = figure();
% Contour plot with curve levels for each point in xseq
[C1, ~] = contour(X, Y, Z, [1, 10, 50, 100, 200, 500, 1000, 1000:1000:9000]);
hold on
% plot of the points in xseq_sd
p = plot(xseq_tn_1(1, :), xseq_tn_1(2, :), '--*');
title('Test 1- Truncated Newton')
legend(p, 'xseq')
hold off;

fig2 = figure();
% Contour plot with curve levels for each point in xseq
[C1, ~] = contour(X, Y, Z, [1, 10, 50, 100, 200, 500, 1000, 1000:1000:9000]);
hold on
% plot of the points in xseq_sd
p = plot(xseq_tn_2(1, :), xseq_tn_2(2, :), '--*');
title('Test 2 - Truncated Newton')
legend(p, 'xseq')
hold off;