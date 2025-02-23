clear
clc
close all

addpath('findiff_general')
addpath('solvers')
addpath('test_builds\exact_builds\')
addpath('test_builds\findiff_builds')

[f1, gradf1, Hessf1] = exact_problem25();


n = 1e5;
x = (rand(n, 1) - 0.5) * 5;

h = sqrt(eps)*1e2;
gradfx = gradf1(x);
gradf_app = zeros(n, 1);
gradf_app(1:2:end) = ...
    0.5*(100 * ( (x(1:2:end) + h).^4 + x(2:2:end).^2 - 2.*x(2:2:end) .* ((x(1:2:end)+h).^2)) - ...
    100 * (x(1:2:end).^2 - x(2:2:end)).^2 + h .* (2*x(1:2:end) - 2 + h) );
gradf_app(2:2:end) = ...
    -50 * h .* (2 .* x(1:2:end).^2 - 2 .* x(2:2:end) - h);
gradf_app = gradf_app / h;

disp(norm(gradf_app - gradfx, "inf"))