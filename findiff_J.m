function Jfx = findiff_J(F, x, h, type, sym, adapt)
%FINDIFF_J
%
% Function that returns an approximation of the jacobian matrix of F
% computed at x
%
% Inputs
% F: R^n -> R^n function
% x: column vector of length n
% h: finite differences "step"
% type: a string of characters. 'fw' -> forward method, 'c' -> centered method
% sym: boolean value. If true, forces Jfx to be symmetric
% adapt: boolean value. If true, when computing the finite differences we use of
%   h*|xbar(i)| as step instead of h
%
% Output
% Jfx: matrix of size n x n, estimation of J_f(x)
%

n = length(x);
Jfx = zeros(n, n);

switch adapt
    case true
        % Adaptative finite difference step
        fstep = @(x, i) h*abs(x(i));
    otherwise
        % Default behaviour
        % Use constant step
        fstep = @(x, i) h;
end

switch type

    case 'Jfw'
        for i = 1:n
            xh = x;
            step = fstep(x, i);
            xh(i) = xh(i) + step;
            Jfx(:, i) = (F(xh) - F(x)) / step;
        end

    case 'Jc'

        for i = 1:n
            xh_plus = x;
            xh_minus = x;
            step = fstep(x, i);
            xh_plus(i) = xh_plus(i) + step;
            xh_minus(i) = xh_minus(i) - step;
            Jfx(:, i) = (F(xh_plus) - F(xh_minus)) / (2 * step);
        end

    otherwise

        for i = 1:n
            xh_plus = x;
            xh_minus = x;
            step = fstep(x, i);
            xh_plus(i) = xh_plus(i) + step;
            xh_minus(i) = xh_minus(i) - step;
            Jfx(:, i) = (F(xh_plus) - F(xh_minus)) / 2 * step;
        end

end

if sym
    Jfx = (Jfx + Jfx') / 2;
end

end

