function [x_stat, g_stat, steps] = bb(f, g, x0, max_iter, tol)
  if nargin() < 3
    error("not enough inputs supplied")
  end
  if not(exist("max_iter"))
    max_iter = 1000
  end
  if not(exist("tol"))
    tol = 10d-10
  end
  l = 1;
  x = [x0];
  d = [-g(x(:,l))];
  alpha = [1];
  while norm(d(:,l)) > tol && l < max_iter
    %% compute next step
    x = [x, x + alpha(l)*d(:,l)];
    %% compute next gradient
    d = [d, -g(x(:,l+1))];
    s = x(:,l+1) - x(:,l);
    y = d(:,l+1) - d(:,l);
    %% compute next alpha
    alpha = [alpha, (s' * s) / (s' * y)];
    l+=1;
    if norm(d(:,l) - d(:,l-1)) < tol
      break
    end
  end
  steps = l;
  x_stat = x;
  g_stat = d;
end
