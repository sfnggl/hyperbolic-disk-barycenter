% https://it.mathworks.com/matlabcentral/answers/787539-steepest-descent-algorithm-in-matlab
clear; clc

% Define the interval constraints
bounds = [3, 9; 3, 9];
% Solve the optimization problem
x_opt = steepest_descent(@obj_func, @grad_func, bounds, 1000, 0.1, 1e-6);

% Print the optimal solution
disp(['Optimal solution: [', num2str(x_opt(1)), ', ', num2str(x_opt(2)), ']']);

% Define the steepest descent method

function x = steepest_descent(func, grad, bounds, max_iter, alpha, tol)
x0 = zeros(size(bounds, 1), 1);
for i = 1:max_iter
    x_prev = x0;
    grad_prev = grad(x_prev);
    x0 = x_prev - alpha * grad_prev;
    % Project onto the feasible set
    for j = 1:length(x0)
        x0(j) = max(bounds(j, 1), min(x0(j), bounds(j, 2)));
    end
    if norm(x0 - x_prev) < tol
        break;
    end
end
x = x0;
end


% Define the objective function to minimize
function f = obj_func(x)
f = x(1).^2 + x(1).*x(2) + 3*x(2).^2;
end

% Define the gradient of the objective function
function g = grad_func(x)
g = [2*x(1)+x(2), x(1)+6*x(2)];
end
