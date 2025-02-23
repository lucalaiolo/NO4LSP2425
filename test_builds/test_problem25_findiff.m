addpath('exact_builds')
addpath('findiff_builds')
rng(1)
adapt = true;
[f, gradf, Hessf] = exact_problem25();
x = (rand(1e5, 1) - 0.5) * 5;
gradfx_approx = findiff_grad_problem25(f, x, sqrt(eps), adapt);
gradfx = gradf(x);
res_grad = abs(gradfx_approx - gradfx);
disp('Problem 25 - Approximation of the Gradient.')
disp(['Adapt = ', num2str(adapt)])
disp(['Maximum error: ', num2str(max(res_grad))])
disp(['Median error: ', num2str(median(res_grad))])
disp(['Mean error: ', num2str(mean(res_grad))])
disp(['Standard deviation: ', num2str(std(res_grad))])

Hessfx_approx = findiff_Hess_problem25(gradf, x, sqrt(eps), adapt);
Hessfx = Hessf(x);

res_diag = abs(full(diag(Hessfx)) - full(diag(Hessfx_approx)));
res_diag_plus = abs(full(diag(Hessfx, 1)) - full(diag(Hessfx_approx, 1)));
% res_diag_minus == res_diag_plus by symmetry
disp('Problem 25 - Approximation of the Hessian.')
disp(['Adapt = ', num2str(adapt)])
disp(['Maximum error on the diagonal: ', num2str(max(res_diag))]);
disp(['Median error: ', num2str(median(res_diag))]);
disp(['Mean error: ', num2str(mean(res_diag))]);
disp(['Standard deviation: ', num2str(std(res_diag))]);
disp(' ')
disp(['Maximum error on the upper diagonal: ', num2str(max(res_diag_plus))]);
disp(['Median error: ', num2str(median(res_diag_plus))]);
disp(['Mean error: ', num2str(mean(res_diag_plus))]);
disp(['Standard deviation: ', num2str(std(res_diag_plus))]);
disp(' ')