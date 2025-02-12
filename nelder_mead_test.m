%% Defining the variables for the test

clear
clc
close all

addpath('solvers\')

% Rosenbrock function
f = @(x) 100 * (x(2,:) - x(1,:).^2).^2 + (1 - x(1,:)).^2;

x0 = [-1.2; 1];

% Default choices for the parameters
rho = 1;
chi = 2;
gamma = 0.5;
sigma = 0.5;
tolFun = 1e-4;
tolX = 1e-4;
MaxIter = 200 * length(x0);
MaxFunEvals = 400 * length(x0);
c = 0.05;

%% Run the tests

[sol, iter, n_fun_evals, f_seq, transformation, best_seq] = ...
    nelder_mead(f, x0, tolFun, tolX, MaxIter, MaxFunEvals, ...
        c, rho, chi, gamma, sigma);

%% Plots

[X, Y] = meshgrid(linspace(-6,6,1000), linspace(-6,6,1000));
f_meshgrid = @(X,Y)reshape(f([X(:),Y(:)]'),size(X));
Z = f_meshgrid(X, Y);

fig1 = figure();
s = 1:10;
[M, ~] = contour(X, Y, Z, [10.^-s 1 10 100 200 500 1000 2000 4000], '-k');
hold on
% plot of the points
plot([x0(1) best_seq(1, :)], [x0(2) best_seq(2,:)], 'r--*')
xlabel('x')
ylabel('y')
title('Nelder-Mead with starting point x0 = [-1.2, 1]');
grid on;
hold off

fig2 = figure();
semilogy([f(x0) f_seq], "bo", "MarkerSize", 3, "MarkerFaceColor", "b");
xlabel('Iterations')
ylabel('Function value')
grid on;
title(['Final function value: ' num2str(f_seq(end))])
hold off;

fig3 = figure();
histogram(categorical(transformation));
title('Transformation counts (x0 = [-1.2, 1])')
hold off;