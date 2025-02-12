function [f, gradf, Hessf] = exact_problem32()
% exact_problem32
%
% Computes the generalized Broyden tridiagonal function in R^n (n > 2), 
% its gradient and its Hessian (in an exact way)
%
% Outputs
% f: function handle that takes column vectors in R^n and returns the
%   function computed at that point
% gradf: function handle that takes column vectors in R^n and returns the
%   gradient of the function computed at that point
% Hessf: function handle that takes column vectors and returns the Hessian
%   matrix of the function computed at that point 
%   (specifically, Hessf(x) is a n x n sparse matrix)
%

fk = @(x) [
    (3 - 2 * x(1,:)) .* x(1,:) + 1 - x(2,:);
    (3 - 2 .* x(2:end-1,:)) .* x(2:end-1,:) + 1 - x(1:end-2,:) - x(3:end,:);
    (3 - 2 * x(end,:)) .* x(end,:) + 1 - x(end-1,:);
];

f = @(x) 0.5 * sum(fk(x).^2)';

gradf = @(x) compute_gradient(x, fk);

% Remark: the generalized Broyden tridiagonal function is smooth
% As a consequence, Hessf is symmetric

Hessf = @(x) compute_Hessian(x, fk);

end

function gradf_build = compute_gradient(x, fk)

fkx = fk(x);

gradf_build = [
    fkx(1) * (3 - 4 * x(1)) - fkx(2);
    - fkx(1:end-2) + fkx(2:end-1) .* (3 - 4 .* x(2:end-1)) - fkx(3:end);
    - fkx(end-1) + (3 - 4 * x(end)) * fkx(end);
];

end

function Hessf_build = compute_Hessian(x, fk)

n = length(x);
fkx = fk(x);

d = [
    (3 - 4 * x(1))^2 - 4 * fkx(1) + 1;
    (3 - 4 * x(2:end-1)).^2 - 4 .* fkx(2:end-1) + 2;
    (3 - 4 * x(end))^2 - 4 * fkx(end) + 1;
];

d_plus = -6 + 4 * x(1:end-1) + 4 * x(2:end); 

d_plus_plus = ones(n - 2, 1);

Hessf_build = spdiags( ...
    [[d_plus_plus;0;0] [d_plus; 0] d [0; d_plus] [0;0;d_plus_plus]], ...
    -2:2, n, n);

end
