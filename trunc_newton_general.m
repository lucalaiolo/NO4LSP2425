function [xk, fk, gradfk_norm, k, xseq, btseq, flag, pcg_iterseq] = ...
    trunc_newton_general(x0, f, gradf, Hessf, kmax, tolgrad, pcg_maxit, ...
        fterms, c1, rho, btmax, FDgrad, FDHess, h, adapt)
% TRUNCATED NEWTON METHOD WITH BACKTRACKING
%
% INPUTS
% x0: column vector of n elements, starting point
% f: function handle to be minimized
% gradf: function handle representing the gradient of f
% Hessf: function handle representing the hessian of f
% kmax: maximum number of iterations allowed
% tolgrad: gradient tolerance 
% pcg_maxit: maximum number of iterations for pcg solver
% fterms: forcing terms for the inexact newton method
% c1: factor c1 for the Armijo condition
% rho: factor less than 1, used to reduce alpha (fixed, for simplicity)
% btmax: maximum number of backtracks
% FDgrad: string that if it is 'fw' or 'c' the newton method will use the
%   corresponding methods to approximate the gradient. Otherwise, the
%   function handle gradf is used
% FDHess: string that if it is 'fw' or 'c' the newton method will use the
%   corresponding methods to approximate the hessian. If it is 'MF', matrix free
%   implementation. Otherwise, the function handle Hessf is used. I
% h: finite differences "step"
% adapt: boolean value. To understand more about this, read findiff_grad
%
% OUTPUTS
% xk: last vector xk computed by the iterative method
% fk: value of function f at xk
% gradfk_norm: norm of the gradient of f computed at xk
% k: index value of the last step executed by the method before stopping
% xseq: matrix n-times-k whose column j is the j-th vector generated by
%   the method
% btseq: row vector of length k that contains the number of backtracks
%   at each step
% flag: extra output variable (string) describing failure or convergence
%   of  the method (and the reasons if there is no convergence)
% pcgiterseq: the number of pcg-iterations executed at each main iteration
%   for computing the direction

switch FDgrad
    case 'fw'
        % OVERWRITE gradf WITH F. HANDLE THAT USES findiff_grad
        % [with option 'fw']
        gradf = @(x) findiff_grad(f, x, h, 'fw', adapt);
    case 'c'
        % OVERWRITE gradf WITH F. HANDLE THAT USES findiff_grad
        % [with option 'c']
        gradf = @(x) findiff_grad(f, x, h, 'c', adapt);
end

% IN THE CASE OF APPROXIMATED GRADIENT, Hessf CAN BE APPROXIMATED ONLY WITH
% findiff_Hess

if isequal(FDgrad, 'fw') || isequal(FDgrad, 'c')
    switch FDHess
        case 'c'
            Hessf = @(x) findiff_Hess(f, x, sqrt(h), adapt);
        case 'MF'
            % MATRIX FREE IMPLEMENTATION
            Hessf_pk = @(x, p) (gradf(x + h*p)-gradf(x)) / h;
    end
else 
    switch FDHess
        case 'c'
            % OVERWRITE Hessf WITH F. HANDLE THAT USES findiff_Hess
            Hessf = @(x) findiff_Hess(f, x, sqrt(h), adapt);
        case 'Jfw'
            % OVERWRITE Hessf WITH F. HANDLE THAT USES findiff_J
            % [with option 'Jfw']
            Hessf = @(x) findiff_J(gradf, x, h, 'Jfw', true, adapt);
        case 'Jc'
            % OVERWRITE Hessf WITH F. HANDLE THAT USES findiff_J
            % [with option 'Jc']
            Hessf = @(x) findiff_J(gradf, x, h, 'Jc', true, adapt);
        case 'MF'
            % MATRIX FREE IMPLEMENTATION
            Hessf_pk = @(x, p) (gradf(x + h*p)-gradf(x)) / h;
    end
end

% Initialization
xk = x0;
fk = f(x0);
gradfk = gradf(xk);

xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);
pcg_iterseq = zeros(1, kmax);

gradfk_norm = norm(gradfk);
k = 0;

farmijo = @(fk, alpha, c1_gradfk_pk) ...
    fk + alpha * c1_gradfk_pk;

while k < kmax && gradfk_norm > tolgrad
    
    % Compute the descent direction as solution of
    % Hessf(xk) p = - gradf(xk)

    switch FDHess
        case 'MF'
            % Initialization
            pk = -gradfk;
            Hpk = Hessf_pk(xk, pk);
            rk = -gradfk - Hpk;
            relres = norm(rk) / gradfk_norm;
            fterms_k = fterms(k, gradfk_norm); 
            j = 0;
            dk = rk;
            
            p_proposed = pk;
            while j < pcg_maxit && relres > fterms_k
                zk = Hessf_pk(xk, dk);
                
                % !!! zk might be close to zero !!!
                if norm(zk) < eps
                    dk = rk;
                    zk = Hessf_pk(xk, dk) + sqrt(eps)*dk;
                end

                alphak_cg = (rk' * dk) / (dk' * zk);
                p_proposed = pk + alphak_cg * dk;
                
                % Check curvature condition
                if p_proposed' * Hessf_pk(xk, p_proposed) <= 0
                    % Exit the loop
                    break;
                end
        
                % Else, accept the proposal
                pk = p_proposed;
        
                rk = rk - alphak_cg * zk;
                betak = -(rk' * zk) / (dk' * zk);
                dk = rk + betak * dk;
            
                relres = norm(rk) / gradfk_norm;
                j = j + 1;
            end
            
            pcg_iterseq(k + 1) = j;

        otherwise
            
            % Initialization
            pk = -gradfk;
            Hessfk = Hessf(xk);
            hessfk_pk = Hessfk * pk;
            rk = -gradfk - hessfk_pk;
            relres = norm(rk) / gradfk_norm;
            fterms_k = fterms(k, gradfk_norm);  
            
            j = 0;
            dk = rk;
            
            p_proposed = pk; 
        
            while j < pcg_maxit && relres > fterms_k
                zk = Hessfk * dk;
                alphak_cg = (rk' * dk) / (dk' * zk);
                p_proposed = pk + alphak_cg * dk;
                
                % Check curvature condition
                if p_proposed' * Hessfk * p_proposed <= 0
                    % Exit the loop
                    break;
                end
        
                % Else, accept the proposal
                pk = p_proposed;
        
                rk = rk - alphak_cg * zk;
                betak = -(rk' * zk) / (dk' * zk);
                dk = rk + betak * dk;
            
                relres = norm(rk) / gradfk_norm;
                j = j + 1;
            end 
            
            pcgiterseq(k+1) = j;

    end

    % Backtracking

    alpha = 1;
    bt = 0;
    
    xnew = xk + alpha * pk;
    fnew = f(xnew);

    c1_gradfk_pk = c1 * gradfk' * pk;

    while bt < btmax && fnew > farmijo(fk, alpha, c1_gradfk_pk)
        alpha = rho * alpha;
        xnew = xk + alpha * pk;
        fnew = f(xnew);
        bt = bt + 1;
    end
    
    if bt == btmax && fnew > farmijo(fk, alpha, c1_gradfk_pk)
        disp('Backtracking strategy could not find a suitable steplength');
        k = k + 1;
        btseq(k) = bt; 
        break
    end

    xk = xnew;
    fk = fnew;
    gradfk = gradf(xk);
    gradfk_norm = norm(gradfk);

    k = k + 1;

    xseq(:, k) = xk;
    btseq(k) = bt;

end

flag = ['Procedure stopped in ', num2str(k), ... 
    ' steps, with gradient norm ', num2str(gradfk_norm)];
xseq = xseq(:, 1:k);
xseq = [x0, xseq];
btseq = btseq(1:k);
pcg_iterseq = pcg_iterseq(1:k);

end