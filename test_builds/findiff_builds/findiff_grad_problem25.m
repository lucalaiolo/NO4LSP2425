function gradfx = findiff_grad_problem25(f, x, h, adapt)
% findiff_grad_problem25
%
% Computes an approximation of the gradient of the extended Rosenbrock function
%
% Inputs:
% f: function handle that takes column vectors in R^n and returns
%   scalar values.
% x: column vector in R^n, point at which we want to evaluate the
%   approximation of the gradient
% h: finite difference "step"
% adapt: boolean value. If true, finite difference steps are computed
%   taking into account abs(x(i))
%
% Output:
% gradfx: approximation of gradf(x)
%
% !!! THIS FUNCTION WORKS ONLY IF n = length(x) IS EVEN !!!

switch adapt
    case true
        step = h * abs(x);
    otherwise
        step = h * ones(length(x), 1);
end
step(step < 1e-12) = h; % avoid dividing by excessively small values when 
% computing the finite differences

n = length(x);
gradfx = zeros(n, 1);
gradfx(1:2:end) = ...
    0.5*(100 * ( (x(1:2:end) + step(1:2:end)).^4 + x(2:2:end).^2 - 2.*x(2:2:end) .* ((x(1:2:end)+step(1:2:end)).^2)) - ...
    100 * (x(1:2:end).^2 - x(2:2:end)).^2 + step(1:2:end) .* (2*x(1:2:end) - 2 + step(1:2:end)) );
gradfx(2:2:end) = ...
    -50 * step(2:2:end) .* (2 .* x(1:2:end).^2 - 2 .* x(2:2:end) - step(2:2:end));
gradfx = gradfx ./ step;

end

