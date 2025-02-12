function Hessfx = findiff_Hess_problem25(gradf, x, h, adapt)
% findiff_Hess_problem25
%
% Computes an approximation of the Hessian of the extended Rosenbrock function
% f whose gradient (or an approximation of it) is gradf.
%
% Inputs:
% gradf: function handle that takes column vectors in R^n and returns
%   column vectors in R^n. In particular, gradf represents the gradient of
%   the function f (be it exact or not, i.e approximated using finite
%   differences)
% x: column vector in R^n, point at which we want to evaluate the
%   approximation of the Hessian
% h: finite difference "step"
% adapt: boolean value. If true, finite difference steps are computed
%   taking into account abs(x(i))
%
% Output:
% Hessfx: approximation of Hessf(x)
%
% !!! THIS FUNCTION WORKS ONLY IF n = length(x) IS EVEN !!!

% With a total of 3 gradf evaluation, we can compute Hessfx


switch adapt
    case true
        step = h * abs(x);
    otherwise
        step = h * ones(length(x), 1);
end
step(step < 1e-12) = h; % avoid dividing by excessively small values when 
% computing the finite differences

n = length(x);
gradfx = gradf(x);

ix_1 = 1:2:n;
ix_2 = 2:2:n;

xh_1 = x;
xh_2 = x;

xh_1(ix_1) = xh_1(ix_1) + step(ix_1);
xh_2(ix_2) = xh_2(ix_2) + step(ix_2); % Perturbed x

main_diag = zeros(n, 1);
lower_diag = zeros(n - 1, 1);
upper_diag = zeros(n - 1, 1); % Diagonals of the Hessian

gradxh_1 = gradf(xh_1);
gradxh_2 = gradf(xh_2);

main_diag(ix_1) = gradxh_1(ix_1);
main_diag(ix_2) = gradxh_2(ix_2);

lower_diag(ix_1) = gradxh_1(ix_2);

upper_diag(ix_1) = gradxh_2(ix_1);

main_diag = (main_diag - gradfx) ./ step;
upper_diag(ix_1) = (upper_diag(ix_1) - gradfx(ix_1)) ./ step(ix_2);
lower_diag(ix_1) = (lower_diag(ix_1) - gradfx(ix_2)) ./ step(ix_1);

Hessfx = spdiags([[lower_diag; 0] main_diag [0; upper_diag]], -1:1, n, n);

Hessfx = (Hessfx + Hessfx') / 2; % Force symmetry

end
