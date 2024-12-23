function [f, gradf, Hessf] = exact_ext_rosenbrock()
% exact_ext_rosenbrock
%
% Computes the extended Rosenbrock function in R^n (n > 2, n EVEN), its gradient 
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

f = @(x) 0.5 * (sum(100 .* (x(1:2:end).^2 - x(2:2:end)).^2) + ... 
        sum((x(1:2:end) - 1).^2));

gradf = @(x) compute_gradient(x);

% Remark: the extended Rosenbrock function is smooth
% As a consequence, Hessf is symmetric

Hessf = @(x) compute_Hessian(x);

end

function gradf_build = compute_gradient(x)

n = length(x);
gradf_build = zeros(n, 1);

gradf_build(1:2:n) = 200 .* x(1:2:n).^3 - 200 .* x(1:2:n) .* x(2:2:n) + x(1:2:n) - 1;
gradf_build(2:2:n) = -100 .* x(1:2:(n-1)).^2 + 100 .* x(2:2:n);

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


