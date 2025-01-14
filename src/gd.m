function [xs, ds, steps] = gd(f,g,x0, max_iter, tol)
  % Gradient Descent with Armijo's Rule for inexact line search
  % It uses a small starting α to offset the rapidly varying derivative
  % Iteratively increases the step taken to update the x with the gradient
  %
  l = 1;
  x = [x0];
  d = [-g(x(:,l))];
  alpha = [armijo(f,g,x(:,l))];
  %% Iterate up to a numerical bound or convergence
  while norm(d(:,l)) > tol && l < max_iter
    % compute gradient
    x_next = x(:,l) + alpha(:,l)*d(:,l);
    if (norm(x_next,1) >= 1)
      break
    end
    x = [x, x_next];
    d = [d,-g(x(:,l+1))];
    % find new alpha via Armijo's method
    % by which chosing α s.t. max[f(x + α*d) < f(x) - σ*d'*d*α]
    % achieves convergence to the result
    if norm(d(:,l+1) - d(:,l)) <= tol
      l = l+1;
      break
    end
    alpha = [alpha,armijo(f,g,x(:,l+1))];
    l = l+1;
  end
  steps = l;
  xs = x;
  ds = d;
end
