function [best, k, n_fun_evals, f_seq, transformation, best_seq] = ...
        nelder_mead(f, x0, tolFun, tolX, MaxIter, MaxFunEvals, ...
            c, rho, chi, gamma, sigma)
%
% Nelder-Mead method
%
% INPUTS
% f: function handle that takes matrices as input and returns vectors. 
%   Specifically, the input matrix should be structured so that its columns 
%   represent the points at which we want to evaluate f.
% x0: starting point, vector in R^n
% tolFun: termination tolerance on the function value, a positive scalar. 
%   Unlike other solvers, fminsearch stops when it satisfies both TolFun 
%   and TolX (in accordance to fminsearch).
% tolX: termination tolerance on x, a positive scalar. Unlike other solvers, 
%   fminsearch stops when it satisfies both TolFun and TolX (in accordance to
%   fminsearch).
% MaxIter: maximum number of iterations allowed, a positive integer
% MaxFunEvals: maximum number of function evaluations allowed, a positive
%   integer. When evaluating matrices of points, the number of function
%   evaluations increases by the number of columns
% c: scalar value. It is used to initialize the simplex. Given x0 as
%   starting point, the starting simplex is simply [x0 x0 + c * eye(n)] where
%   n is equal to length(x0).
% rho: reflection parameter, scalar positive value
% chi: expansion parameter, scalar positive value greater than 1
% gamma: contraction parameter, 0 < gamma < 1
% sigma: shrinkage parameter, 0 < sigma < 1
%
% OUTPUTS
% best: found solution, column vector
% k: index value of the last step executed by the method before stopping
% n_fun_evals: total number of function evaluations 
% f_seq: vector in R^{k+1}, f_seq(i) is the function value at the i-th step
% transformation: cell array containing the sequence of transformations
%   applied to the initial simplex
%
n = length(x0);
simplex = [x0 x0 + c * eye(n)]; % standard initialization
k = 0;
f_seq = zeros(1, MaxIter);
transformation = cell(1, MaxIter);
best_seq = zeros(n, MaxIter); % can be used if n is not too large

fx = f(simplex);
[fx, order] = sort(fx); % f must be able to take matrices as input
simplex = simplex(:, order); % ordered initial simplex
fx = fx(order);
res_x = max(vecnorm(simplex(:, 1) - simplex(:, 2:end), "inf"));
res_fun = abs(fx(1) - fx(end));
n_fun_evals = n + 1;

while ~(k >= MaxIter || (res_x < tolX && res_fun < tolFun) || n_fun_evals > MaxFunEvals)
    % Sort simplex based on function values
    
    best = simplex(:, 1);
    worst = simplex(:, end);

    % Compute centroid excluding worst point
    x_bar = mean(simplex(:, 1:end-1), 2);

    % Reflection step
    x_r = x_bar + rho * (x_bar - worst);
    f_reflection = f(x_r);
    n_fun_evals = n_fun_evals + 1;
    f_worst = fx(end);
    if f_reflection < fx(1)
        % Expansion step
        x_e = x_bar + chi * (x_r - x_bar);
        f_expansion = f(x_e);
        n_fun_evals = n_fun_evals + 1;
        if f_expansion < f_reflection
            % x_e is the new best point
            simplex(:, 2:end) = simplex(:, 1:end-1);
            simplex(:, 1) = x_e;
            fx(2:end) = fx(1:end-1);
            fx(1) = f_expansion;
            transformation{k+1} = 'expansion';
        else
            % x_r is the new best point
            simplex(:, 2:end) = simplex(:, 1:end-1);
            simplex(:, 1) = x_r;
            fx(2:end) = fx(1:end-1);
            fx(1) = f_reflection;
            transformation{k+1} = 'reflection';
        end
    elseif f_reflection >= fx(end-1) % f(second_worst)
        % Contraction step
        if f_reflection < f_worst
            x_c = x_bar + gamma * ( x_r - x_bar );
        else
            x_c = x_bar + gamma * ( worst - x_bar );
        end
        f_contraction = f(x_c);
        n_fun_evals = n_fun_evals + 1;
        if f_contraction < f_worst
            % Accept
            simplex(:, end) = x_c;
            fx(end) = f_contraction;
            [fx, order] = sort(fx);
            simplex = simplex(:, order);
            transformation{k+1} = 'contraction';
        else
            % Shrink step
            transformation{k+1} = 'shrinkage';
            simplex(:, 2:end) = best + sigma * (simplex(:, 2:end) - best);
            fx(2:end) = f(simplex(:, 2:end));
            [fx, order] = sort(fx);
            simplex = simplex(:, order);
            n_fun_evals = n_fun_evals + n;
        end
    else
        % Accept x_r
        simplex(:, end) = x_r;
        fx(end) = f_reflection;
        [fx, order] = sort(fx);
        simplex = simplex(:, order);
        transformation{k+1} = 'expansion';
    end

    % Update res_x, res_fun and f_seq
    f_seq(k + 1) = fx(1);
    res_x = max(vecnorm(simplex(:, 1) - simplex(:, 2:end), "inf"));
    res_fun = abs(fx(1) - fx(end));
    best_seq(:, k + 1) = simplex(:, 1);
    k = k + 1;
end

best = simplex(:, 1);
f_seq = f_seq(1:k);
transformation = transformation(1:k);
best_seq = best_seq(:, 1:k);
end