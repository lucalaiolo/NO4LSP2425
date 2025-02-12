function [pk, k, relres] = ...
    pcg_curvature_condition(A, b, kmax, tol, M1, M2)
%
% Preconditioned Conjugate Gradient Method + negative curvature condition 
% check for the truncated Newton method
%
% INPUTS:
% A : function handle. A(x) computes the matrix-vector product A * x
% b : column vector of n elements, known terms of the system (-gradfk)
% kmax : positive integer, maximum number of steps
% tol : positive real value, tolerance for relative residual
% M1, M2 : matrices such that M = M1 * M2 is the chosen preconditioner
%   (typically, M1 = L where L = ichol(A) and M2 = L')
%
% OUTPUTS:
% pk : column vector of n elements, solution computed
% k : integer, number of iterations taken to compute xk
% relres : relative residual
%

zk = zeros(length(b), 1);
zk_old = b;
rk = b;
yk = M2 \ (M1 \ rk);
dk = yk;

k = 0;

norm_b = norm(b);
relres = norm(rk) / norm_b;

while k < kmax && relres > tol
    gk = A(dk);
    % Curvature condition check
    if dk' * gk <= 0
        if k==0 || k == 1
            % Exit with the steepest descent direction
            pk = b;
            return;
        else
            % Exit with the last valid iterate
            pk = zk_old;
            return;
        end
    end
    alphak = (rk' * yk) / (dk' * gk);
    zk_old = zk;
    zk = zk + alphak * dk;
    
    % Update beta
    denom = rk' * yk;
    rk = rk - alphak * gk;
    yk = M2 \ (M1 \ rk);
    num = rk' * yk;
    betak = num / denom;

    dk = yk + betak * dk;

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
