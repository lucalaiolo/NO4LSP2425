function Hessfx = findiff_Hess_tridiagonal(gradf, x, h)
% findiff_Hess_tridiagonal
%
% Computes an approximation of the Hessian of a function f whose gradient
% (or an approximation of it) is gradf and whose Hessian is tridiagonal
% (with diagonals on positions -1, 0, 1)
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


% With a total of 4 gradf evaluation, we can compute Hessfx

n = length(x);
gradfx = gradf(x);

ix_1 = 1:3:n;
ix_2 = 2:3:n;
ix_3 = 3:3:n;

xh_1 = x;
xh_2 = x;
xh_3 = x;

xh_1(ix_1) = xh_1(ix_1) + h;
xh_2(ix_2) = xh_2(ix_2) + h;
xh_3(ix_3) = xh_3(ix_3) + h; % Perturbed x

main_diag = zeros(n, 1);
lower_diag = zeros(n - 1, 1);
upper_diag = zeros(n - 1, 1); % Diagonals of the Hessian

gradxh_1 = gradf(xh_1);
gradxh_2 = gradf(xh_2);
gradxh_3 = gradf(xh_3);

main_diag(ix_1) = gradxh_1(ix_1);
main_diag(ix_2) = gradxh_2(ix_2);
main_diag(ix_3) = gradxh_3(ix_3);

lower_diag(1:3:min(n - 1, ix_1(end))) = gradxh_1(2:3:n);
lower_diag(2:3:min(n - 1, ix_2(end))) = gradxh_2(3:3:n);
lower_diag(3:3:min(n - 1, ix_3(end))) = gradxh_3(4:3:n);

upper_diag(3:3:(n-1)) = gradxh_1(3:3:(ix_1(end) - 1));
upper_diag(1:3:(n-1)) = gradxh_2(1:3:(ix_2(end) - 1));
upper_diag(2:3:(n-1)) = gradxh_3(2:3:(ix_3(end) - 1));

main_diag = (main_diag - gradfx) / h;
upper_diag = (upper_diag - gradfx(1:n-1)) / h;
lower_diag = (lower_diag - gradfx(2:n)) / h;

Hessfx = spdiags([[lower_diag; 0] main_diag [0; upper_diag]], -1:1, n, n);

Hessfx = (Hessfx + Hessfx') / 2; % Force symmetry

end
