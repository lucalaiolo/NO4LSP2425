addpath('exact_builds')
addpath('findiff_builds')
rng(1)
adapt = true;
[f, gradf, Hessf] = exact_problem31();
x = (rand(1e5, 1) - 0.5) * 5;
Hessfx_approx = findiff_Hess_problem31_32(gradf, x, sqrt(eps), adapt);
Hessfx = Hessf(x);

res_diag = abs(full(diag(Hessfx)) - full(diag(Hessfx_approx)));
res_diag_plus = abs(full(diag(Hessfx, 1)) - full(diag(Hessfx_approx, 1)));
res_diag_plus_plus = abs(full(diag(Hessfx, 2)) - full(diag(Hessfx_approx, 2)));
% res_diag_minus == res_diag_plus by symmetry
% res_diag_minus_minus == res_diag_plus_plus by symmetry
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
disp(['Maximum error on the upper-upper diagonal: ', num2str(max(res_diag_plus_plus))]);
disp(['Median error: ', num2str(median(res_diag_plus_plus))]);
disp(['Mean error: ', num2str(mean(res_diag_plus_plus))]);
disp(['Standard deviation: ', num2str(std(res_diag_plus_plus))]);