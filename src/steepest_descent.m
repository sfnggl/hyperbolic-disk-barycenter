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
  while l < max_iter
    %% compute gradient
    x = [x, x(:,l) + alpha(:,l)*d(:,l)];
    d = [d,-g(x(:,l+1))];
    alpha = [alpha,armijo(f,g,x(:,l+1))];
    l += 1;
    if norm(d(:,l) - d(:,l-1)) < tol
      break
    end
  end
  l
  x_stat = x;
  g_stat = d;
end
