function [f, gradf, Hessf] = exact_problem25()
% exact_problem25
%
% Computes the extended Rosenbrock function in R^n, its gradient 
% and its Hessian (in an exact way)
%
% Outputs
% f: function handle that takes column vectors in R^n and returns the
%   extended Rosenbrock function computed at that point
% gradf: function handle that takes column vectors in R^n and returns the
%   gradient of the extended Rosenbrock function computed at that point
% Hessf: function handle that takes column vectors and returns the Hessian
%   matrix of the extended Rosenbrock function computed at that point 
%   (specifically, Hessf(x) is a n x n sparse matrix)
%

f = @(x) 0.5 * (sum(100 .* (x(1:2:end-1, :).^2 - x(2:2:end, :)).^2) + ... 
        sum((x(1:2:end-1, :) - 1).^2)) + 0.5 * mod(length(x(:,1)), 2) * 100 * x(end)^4;

gradf = @(x) compute_gradient(x);

% Remark: the extended Rosenbrock function is smooth
% As a consequence, Hessf is symmetric

Hessf = @(x) compute_Hessian(x);

end

function gradf_build = compute_gradient(x)

n = length(x);
p = mod(n, 2);

gradf_build = zeros(n, 1);
switch p
    case 0
        % n is even
        gradf_build(1:2:n) = 200 .* x(1:2:n).^3 - 200 .* x(1:2:n) .* x(2:2:n) + x(1:2:n) - 1;
        gradf_build(2:2:n) = -100 .* x(1:2:(n-1)).^2 + 100 .* x(2:2:n);
    case 1
        % n is odd --> n+1 is even
        gradf_build = [gradf_build; 0];
        x_aug = [x; 0];
        gradf_build(1:2:end) = 200 .* x_aug(1:2:end).^3 - 200 .* x_aug(1:2:end) ...
            .* x_aug(2:2:end) + x_aug(1:2:end) - 1;
        gradf_build(2:2:end) = -100 .* x_aug(1:2:(end-1)).^2 + 100 .* x_aug(2:2:end);
        gradf_build = gradf_build(1:end-1);

end

end

function Hessf_build = compute_Hessian(x)

n = length(x);

d = zeros(n, 1);
d_plus = zeros(n - 1, 1);

d(1:2:n) = 600 .* x(1:2:n).^2 - 200 .* x(2:2:n) + 1;
d(2:2:n) = 100; % Main diagonal

d_plus(1:2:(n-1)) = - 200 .* x(1:2:(n-1)); % Upper diagonal

Hessf_build = spdiags([[d_plus; 0] d [0; d_plus]], -1:1, n, n);

end


