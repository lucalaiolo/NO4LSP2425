function Hessfx = findiff_Hess_diagonal(gradf, x, h)
% findiff_Hess_tridiagonal
%
% Computes an approximation of the Hessian of a function f whose gradient
% (or an approximation of it) is gradf and whose Hessian is diagonal
% (only the main diagonal is non zero)
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


% With a total of 2 gradf evaluation, we can compute Hessfx

n = length(x);
gradfx = gradf(x);

xh = x + h;

gradxh = gradf(xh);

main_diag = (gradxh - gradfx) / h;

Hessfx = sparse(1:n, 1:n, main_diag);

end
