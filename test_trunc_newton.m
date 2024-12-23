%% Defining the variables for the test

clear
clc
close all

seed = min([9876, 1234, 2345]); % to do

d = 3:5;
dim = 10.^d;

for i = 1:length(dim)
    n = dim(i);

    % Rosenbrock suggested starting point
    x0_rosenbrock = ones(n, 1);
    x0_rosenbrock(1:2:end) = -1.2;

    % Extended Rosenbrock suggested starting point
    x0_ext_rosenbrock = x0_rosenbrock; % same value

    % Banded trigonometric problem suggested starting point
    x0_banded = ones(n, 1);

    % Perturbed initial conditions
    x0_rosenbrock_rand = x0_rosenbrock + 2 .* rand(n, 10) - 1;
    x0_ext_rosenbrock_rand = x0_ext_rosenbrock + 2 .* rand(n, 10) - 1;
    x0_banded = x0_banded + 2 .* rand(n, 10) - 1;

    % Solve the problems

    % to do






end