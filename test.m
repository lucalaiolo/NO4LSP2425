[f, gradf, Hessf] = exact_rosenbrock();

x = rand(10000,1);
Hessfx = Hessf(x);

Hessfx_approx = findiff_Hess_tridiagonal(gradf, x, sqrt(eps));
disp(['Error (norm infinity) tridiag: ', num2str(norm(Hessfx - Hessfx_approx, "inf"))]);

[f, gradf, Hessf] = exact_banded_trigon();
x = rand(10000, 1);
Hessfx = Hessf(x);
Hessfx_approx = findiff_Hess_diagonal(gradf, x, sqrt(eps));
disp(['Error (norm infinity) diag: ', num2str(norm(Hessfx - Hessfx_approx, "inf"))]);

[f, gradf, Hessf] = exact_ext_rosenbrock();
x = rand(10000, 1);
Hessfx = Hessf(x);
Hessfx_approx = findiff_Hess_ext_rosenbrock(gradf, x, sqrt(eps));
disp(['Error (norm infinity) extended: ', num2str(norm(Hessfx - Hessfx_approx, "inf"))]);
