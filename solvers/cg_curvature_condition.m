function [pk, k, relres] = cg_curvature_condition(A, b, kmax, tol)
%
% Conjugate Gradient Method + negative curvature condition check for the
% truncated Newton method
%
% INPUTS:
% A : function handle. A(x) computes the matrix-vector product A * x
% b : column vector of n elements, known terms of the system (-gradfk)
% kmax : positive integer, maximum number of steps
% tol : positive real value, tolerance for relative residual
%
% OUTPUTS:
% pk : column vector of n elements, solution computed
% k : integer, number of iterations taken to compute xk
% relres : relative residual
%


zk = zeros(length(b),1); % canon choice for this method
zk_old = b;

rk = b; % = -gradfk
dk = rk;

k = 0;

norm_b = norm(b);
relres = norm(rk) / norm_b;
while k < kmax
    
    yk = A(dk);
    % Curvature condition check
    if dk' * yk <= 0
        if k == 0 || k == 1
            % Exit with the steepest descent direction
            pk = b;
            return;
        else
            % Exit with the last valid iterate
            pk = zk_old;
            return;
        end
    end

    alphak = (rk' * dk) / (dk' * yk);
    zk_old = zk;
    zk = zk + alphak * dk;

    rk = rk - alphak * yk;
    betak = -(rk' * yk) / (dk' * yk);
    dk = rk + betak * dk;

    relres = norm(rk) / norm_b;
    if relres <= tol
        % Exit with the current solution
        pk = zk;
        return;
    end
    k = k + 1;
end
pk = zk;
end
