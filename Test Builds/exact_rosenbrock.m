function [f, gradf, Hessf] = exact_rosenbrock()
% exact_rosenbrock
%
% Computes the Rosenbrock function in R^n (n > 2), its gradient and its
% Hessian (in an exact way)
%
% Outputs
% f: function handle that takes column vectors in R^n and returns the
%   Rosenbrock function computed at that point
% gradf: function handle that takes column vectors in R^n and returns the
%   gradient of the Rosenbrock function computed at that point
% Hessf: function handle that takes column vectors and returns the Hessian
%   matrix of the Rosenbrock function computed at that point (specifically,
%   Hessf(x) is a n x n sparse matrix)
%

f = @(x) (sum(100 .* (x(2:end,:) - x(1:end-1,:).^2).^2 + (1 - x(1:end-1,:)).^2))';
% This way, f can take multiple vectors

gradf = @(x) [
    x(1) * 400 * (x(1)^2 - x(2)) + 2*x(1) - 2;
    x(2:end-1).*(202 + 400 .* (x(2:end-1).^2 - x(3:end))) - 2 - 200 .* x(1:end-2).^2 ;
    200 * (x(end) - x(end-1)^2);
];

%{
202 .* x(2:end-1) - 200 .* x(1:end-2).^2 - 400 .* x(2:end-1) .* x(3:end) - ...
        2 + 400 * x(2:end-1).^3
%}

% Remark: the Rosenbrock function is smooth
% As a consequence, Hessf is symmetric

d = @(x) [
    -400 * x(2) + 1200 * x(1)^2 + 2;
    202 - 400 .* x(3:end) + 1200 .* x(2:end-1).^2;
    200;
]; % This is the diagonal of Hessf

d_plus = @(x) -400 .* x(1:end-1); % Upper diagonal of Hessf
% Remark: this is a column vector of length n - 1

% Since Hessf is symmetric and tridiagonal, this is enough to compute it
Hessf = @(x) spdiags([ [d_plus(x); 0] d(x) [0; d_plus(x)] ], -1:1, length(x), length(x));


end

