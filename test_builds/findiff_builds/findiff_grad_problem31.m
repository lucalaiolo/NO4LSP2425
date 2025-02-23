function gradfx = findiff_grad_problem31(f, x, h, adapt)
% findiff_grad_problem25
%
% Computes an approximation of the gradient of the Broyden tridiagonal function
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

gradfx(1) = 0.5 .* ( ((3 - 2 * (x(1) + step(1))) * (x(1) + step(1)) + 1 - 2 * x(2))^2 - ...
    ((3 - 2 * x(1)) * x(1) + 1 - 2 * x(2))^2 - ...
    step(1) * (2 * x(2) * (3 - 2 * x(2)) + 2 - 2 * x(1) - 4 * x(3) - step(1)));

gradfx(2) = 0.5 .* (-2.*step(2) .* ( 2.* x(1) .* (3 - 2 .* x(1)) + ...
    2 - 4.*x(2) - 2.*step(2)) + ...
    ((3 - 2.* (x(2) + step(2))).*(x(2) + step(2)) + 1 - x(1) - 2.*x(3)).^2 - ...
    ((3 - 2.*x(2)) .* x(2) + 1 - x(1) - 2.*x(3)).^2 - ...
    step(2) .* (2.*x(3).*(3 - 2.*x(3)) + 2 - 2*x(2) - 4.*x(4) - step(2)));

gradfx(3:end-2) = 0.5 .* (-2.*step(3:end-2) .* ( 2.* x(2:end-3) .* (3 - 2 .* x(2:end-3)) + ...
    2 - 2 .* x(1:end-4) - 4.*x(3:end-2) - 2.*step(3:end-2)) + ...
    ((3 - 2.* (x(3:end-2) + step(3:end-2))).*(x(3:end-2) + step(3:end-2)) + 1 - x(2:end-3) - 2.*x(4:end-1)).^2 - ...
    ((3 - 2.*x(3:end-2)) .* x(3:end-2) + 1 - x(2:end-3) - 2.*x(4:end-1)).^2 - ...
    step(3:end-2) .* (2.*x(4:end-1).*(3 - 2.*x(4:end-1)) + 2 - 2.*x(3:end-2) - 4.*x(5:end) - step(3:end-2)));

gradfx(end-1) = 0.5 .* (-2.*step(end-1) .* ( 2.* x(end-2) .* (3 - 2 .* x(end-2)) + ...
    2 - 2*x(end-3) - 4.*x(end-1) - 2.*step(end-1)) + ...
    ((3 - 2.* (x(end-1) + step(end-1))).*(x(end-1) + step(end-1)) + 1 - x(end-2) - 2.*x(end)).^2 - ...
    ((3 - 2.*x(end-1)) .* x(end-1) + 1 - x(end-2) - 2.*x(end)).^2 - ...
    step(end-1) .* (2.*x(end).*(3 - 2.*x(end)) + 2 - 2*x(end-1) - step(end-1)));

gradfx(end) = 0.5 * (-2 * step(end) * ( 2*x(end-1) * (3 - 2*x(end-1)) + 2 - 2*x(end-2) - 4*x(end) - 2*step(end)) + ...
    ((3 - 2.* (x(end) + step(end))).*(x(end) + step(end)) + 1 - x(end-1))^2 - ...
    ((3 - 2.*x(end)) .* x(end) + 1 - x(end-1)).^2);

gradfx = gradfx ./ step;



end

