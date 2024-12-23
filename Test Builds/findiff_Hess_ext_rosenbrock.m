function Hessfx = findiff_Hess_ext_rosenbrock(gradf, x, h)
% findiff_Hess_ext_rosenbrock
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
%
% Output:
% Hessfx: approximation of Hessf(x)
%
% !!! THIS FUNCTION WORKS ONLY IF n = length(x) IS EVEN !!!

% With a total of 3 gradf evaluation, we can compute Hessfx

n = length(x);
gradfx = gradf(x);

ix_1 = 1:2:n;
ix_2 = 2:2:n;

xh_1 = x;
xh_2 = x;

xh_1(ix_1) = xh_1(ix_1) + h;
xh_2(ix_2) = xh_2(ix_2) + h; % Perturbed x

main_diag = zeros(n, 1);
lower_diag = zeros(n - 1, 1);
upper_diag = zeros(n - 1, 1); % Diagonals of the Hessian

gradxh_1 = gradf(xh_1);
gradxh_2 = gradf(xh_2);

main_diag(ix_1) = gradxh_1(ix_1);
main_diag(ix_2) = gradxh_2(ix_2);

lower_diag(ix_1) = gradxh_1(ix_2);

upper_diag(ix_1) = gradxh_2(ix_1);

main_diag = (main_diag - gradfx) / h;
upper_diag(ix_1) = (upper_diag(ix_1) - gradfx(ix_1)) / h;
lower_diag(ix_1) = (lower_diag(ix_1) - gradfx(ix_2)) / h;

Hessfx = spdiags([[lower_diag; 0] main_diag [0; upper_diag]], -1:1, n, n);

Hessfx = (Hessfx + Hessfx') / 2; % Force symmetry

end
