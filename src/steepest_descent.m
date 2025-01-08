function [x_stat, g_stat] = steepest_descent(f, g, x0, max_iter, tol)
  if not(exist("max_iter"))
    max_iter = 1000
  end
  if not(exist("tol"))
    tol = 10d-10
  end
  l = 1;
  x = [x0];
  d = [-g(x(:,l))];
  alpha = [armijo(f,g,x(:,l))];
  %% Iterate up to a numerical bound or convergence
  while norm(d(:,l)) > tol && l < max_iter
    %% compute gradient
    x = [x, x(:,l) + alpha(:,l)*d(:,l)];
    d = [d,-g(x(:,l+1))];
    %% random step to prevent stagnatio
    if ~mod(l,25)
      alpha = [alpha, 0.1];
    else
      alpha = [alpha,armijo(f,g,x(:,l+1))];
    end
    l += 1;
  end
  l
  x_stat = x;
  g_stat = d;
end
