function Hessfx = findiff_Hess(f, x, h, adapt)
%FINDIFF_HESS
% Function that computes and approximation of the hessian of f computed
% at x.
%
% !!! GENERAL PURPOSE FUNCTION. CANNOT BE USED WHEN THE DIMENSION IS LARGE
% IT DOES NOT TAKE INTO ACCOUNT SPARSITY OF THE HESSIAN OR OTHER
% INFORMATION
%
% Inputs:
% f: function handle
% x: point at which an approximation of the hessian of f will be computed
% h: finite difference "step"
% adapt: boolean value. If true, when computing the finite differences we use of
%   h*|xbar(i)| as step instead of h
%
% Output:
% Hessfx: approximation of Hessf(x)
%

n = length(x);
Hessfx = zeros(n, n);

switch adapt
    case true
        % Adaptative finite difference step
        fstep = @(x, i) h*abs(x(i));
    otherwise
        % Default behaviour
        % Use constant step
        fstep = @(x, i) h;
end

for i=1:n
    xh_plus = x;
    xh_minus = x;
    step_i = fstep(x, i);
    if step_i == 0
        step_i = h;
    end
    xh_plus(i) = xh_plus(i) + step_i;
    xh_minus(i) = xh_minus(i) - step_i;

    Hessfx(i,i) = (f(xh_plus) - 2*f(x) + f(xh_minus))/(step_i)^2;
    for j=(i+1):n
        step_j = fstep(x, j);
        if step_j == 0
            step_j = h;
        end
        xh_plus_ij = x;
        xh_plus_ij([i,j]) = xh_plus_ij([i,j]) + [step_i; step_j];
        xh_plus_iterate = x;
        xh_plus_iterate(j) = xh_plus_iterate(j) + step_j;

        Hessfx(i,j) = (f(xh_plus_ij) - f(xh_plus) - f(xh_plus_iterate) + f(x))/ ...
            (step_i * step_j); % In case of adapt = false, the denominator
                               % is equal to h^2
        Hessfx(j,i) = Hessfx(i,j);
    end
end

