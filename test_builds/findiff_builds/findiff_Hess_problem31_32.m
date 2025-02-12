function Hessfx = findiff_Hess_problem31_32(gradf, x, h, adapt)
% findiff_Hess_problem31_32
%
% Computes an approximation of the Hessian of the (extended) generalized
% Broyden tridiagonal function f whose gradient (or an approximation of it)
% is gradf.
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
% !!! THIS FUNCTION WORKS ONLY IF n = length(x) > 4 !!!
%
% With a total of 6 gradf evaluations we can compute Hessfx
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
gradfx = gradf(x);

main_diag = zeros(n, 1);
upper_diag = zeros(n - 1, 1);
upper_upper_diag = zeros(n - 2, 1);
lower_diag = zeros(n - 1, 1);
lower_lower_diag = zeros(n - 2, 1);

% The following operations could have been done with a for cycle. However,
% it is far more readable like that
xh = x;
xh(1:5:end) = xh(1:5:end) + step(1:5:end);
gradxh_1 = gradf(xh);

xh = x;
xh(2:5:end) = xh(2:5:end) + step(2:5:end);
gradxh_2 = gradf(xh);

xh = x;
xh(3:5:end) = xh(3:5:end) + step(3:5:end);
gradxh_3 = gradf(xh);

xh = x;
xh(4:5:end) = xh(4:5:end) + step(4:5:end);
gradxh_4 = gradf(xh);

xh = x;
xh(5:5:end) = xh(5:5:end) + step(5:5:end);
gradxh_5 = gradf(xh);

main_diag(1:5:end) = gradxh_1(1:5:end) - gradfx(1:5:end);
main_diag(2:5:end) = gradxh_2(2:5:end) - gradfx(2:5:end);
main_diag(3:5:end) = gradxh_3(3:5:end) - gradfx(3:5:end);
main_diag(4:5:end) = gradxh_4(4:5:end) - gradfx(4:5:end);
main_diag(5:5:end) = gradxh_5(5:5:end) - gradfx(5:5:end);

main_diag = main_diag ./ step;

upper_diag(1:5:end) = gradxh_2(1:5:(end-1)) - gradfx(1:5:(end-1));
upper_diag(2:5:end) = gradxh_3(2:5:(end-1)) - gradfx(2:5:(end-1));
upper_diag(3:5:end) = gradxh_4(3:5:(end-1)) - gradfx(3:5:(end-1));
upper_diag(4:5:end) = gradxh_5(4:5:(end-1)) - gradfx(4:5:(end-1));
upper_diag(5:5:end) = gradxh_1(5:5:(end-1)) - gradfx(5:5:(end-1));

upper_diag = upper_diag ./ step(2:end);

upper_upper_diag(1:5:end) = gradxh_3(1:5:(end-2)) - gradfx(1:5:(end-2));
upper_upper_diag(2:5:end) = gradxh_4(2:5:(end-2)) - gradfx(2:5:(end-2));
upper_upper_diag(3:5:end) = gradxh_5(3:5:(end-2)) - gradfx(3:5:(end-2));
upper_upper_diag(4:5:end) = gradxh_1(4:5:(end-2)) - gradfx(4:5:(end-2));
upper_upper_diag(5:5:end) = gradxh_2(5:5:(end-2)) - gradfx(5:5:(end-2));

upper_upper_diag = upper_upper_diag ./ step(3:end);

lower_diag(1:5:end) = gradxh_1(2:5:end) - gradfx(2:5:end);
lower_diag(2:5:end) = gradxh_2(3:5:end) - gradfx(3:5:end);
lower_diag(3:5:end) = gradxh_3(4:5:end) - gradfx(4:5:end);
lower_diag(4:5:end) = gradxh_4(5:5:end) - gradfx(5:5:end);
lower_diag(5:5:end) = gradxh_5(6:5:end) - gradfx(6:5:end);

lower_diag = lower_diag ./ step(1:end-1);

lower_lower_diag(1:5:end) = gradxh_1(3:5:end) - gradfx(3:5:end);
lower_lower_diag(2:5:end) = gradxh_2(4:5:end) - gradfx(4:5:end);
lower_lower_diag(3:5:end) = gradxh_3(5:5:end) - gradfx(5:5:end);
lower_lower_diag(4:5:end) = gradxh_4(6:5:end) - gradfx(6:5:end);
lower_lower_diag(5:5:end) = gradxh_5(7:5:end) - gradfx(7:5:end);

lower_lower_diag = lower_lower_diag ./ step(1:end-2);

Hessfx = spdiags([ ...
        [lower_lower_diag;0;0] ...
        [lower_diag;0] ...
        main_diag ...
        [0;upper_diag] ...
        [0;0;upper_upper_diag] ...
    ], -2:2, n, n);

Hessfx = (Hessfx + Hessfx') ./ 2; % Force symmetry

end

