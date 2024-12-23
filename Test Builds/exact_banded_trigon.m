function [f, gradf, Hessf] = exact_banded_trigon()
% exact_banded_trigon
%
% Computes the function of the banded trigonometric problem in R^n (n > 2),
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

f = @(x) 1 - cos(x(1)) - sin(x(2)) + sum((2:(n-1))' * (1 - cos(x(2:(n-1))) + ...
        sin(x(1:(n-2))) - sin(x(3:n)) )) + n * (1 - cos(x(n))  + sin(x(n-1))); 

gradf = @(x) [
    (1:(length(x)-1))' .* sin(x(1:(end-1))) + 2 .* cos(x(1:(end-1)));
    cos(x(end)) + length(x) * sqrt(2) * sin(x(end) - pi/4);
];

% The Hessian is diagonal

main_diag = @(x) [
    (1:(length(x) - 1))' .* cos(x(1:(end-1))) - 2 .* sin(x(1:(end-1)));
    -sin(x(end)) + length(x) * sqrt(2) * cos(x(end) - pi/4);
];

Hessf = @(x) sparse(1:length(x), 1:length(x), main_diag(x));

end
