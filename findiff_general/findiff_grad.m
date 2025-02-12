function gradfx = findiff_grad(f, x, h, type, adapt)
% FINDIFF_GRAD 
% Function that computes an approximation of the gradient of f at point x
%
% Inputs
% f: function handle
% x: point at which we want to approximate gradf
% h: the finite difference "step"
% type: a string of characters. 'fw' -> forward method, 'c' -> centered method
% adapt: boolean value. If true, when computing the approximation of
%   df/dx_i(xbar) we use h*|xbar(i)| instead of h
%
% Output
% gradfx: column vector, estimation of gradf(x)
%

switch adapt
    case true
        % Adaptative finite difference step
        fstep = @(x, i) h*abs(x(i));
    otherwise
        % Default behaviour
        % Use constant step
        fstep = @(x, i) h;
end

n = length(x);
gradfx = zeros(n, 1);

switch type
    case 'fw'
        % Forward differences approximation
        fx = f(x);
        for i = 1:n
            xh = x;
            step = fstep(x, i); % we should check whether x(i) is zero or not
            % x(i) is zero if and only if step is zero (if adapt is true)
            if step == 0
                step = h; % we must avoid dividing by zero
            end
            xh(i) = xh(i) + step;
            gradfx(i) = (f(xh) - fx) / step;
        end

    case 'c'
        % Centered differences approximation
        for i = 1:n
            xh_plus = x;
            xh_minus = x;
            step = fstep(x, i); % could have used xh_minus, equivalent
            if step == 0
                step = h;
            end
            xh_plus(i) = xh_plus(i) + step;
            xh_minus(i) = xh_minus(i) - step;
            gradfx(i) = (f(xh_plus) - f(xh_minus)) / (2*step);
        end

    otherwise
        % Default option
        % Forward differences approximation
        fx = f(x);
        for i = 1 : n
            xh = x;
            step = fstep(x, i);
            if step == 0
                step = h;
            end
            xh(i) = xh(i) + step;
            gradfx(i) = (f(xh) - fx) / step;
        end

end

end

