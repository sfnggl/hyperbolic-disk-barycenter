function [xs, ds, steps] = wolfe(f,g,x0, max_iter, tol)
  % Steepest Descent with modified Armijo's Rule for inexact line search
  %
  % Iteratively reduces the step taken to update the x with the gradient
  %
  l = 1;
  x = [x0];
  d = [-g(x(:,l))];
  alpha = [modiarmijo(f,g,x(:,l))];
  %% Iterate up to a numerical bound or convergence
  while norm(d(:,l)) > tol && l < max_iter
    x_next = x(:,l) + alpha(:,l)*d(:,l);
    if (norm(x_next,1) >= 1)
      break
    end
    x = [x, x_next];
    d = [d,-g(x(:,l+1))];
    if norm(d(:,l+1) - d(:,l)) <= tol
      l = l+1;
      break
    end
    % find new alpha via modified Armijo's condition
    % such that α
    % achieves convergence to the result
    alpha = [alpha,modiarmijo(f,g,x(:,l+1))];
    l = l+1;
  end
  steps = l;
  xs = x;
  ds = d;
end
